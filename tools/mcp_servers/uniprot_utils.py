import requests
import time
from tqdm import tqdm
import sys
import multiprocessing
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cache_manager
from well_studied_kinase_utils import get_well_studied_human_kinase_ids, get_domain_seq_from_id

#from cms_utils import cache_manager
#from cms_utils.well_studied_kinase_utils import get_well_studied_human_kinase_ids, get_domain_seq_from_id

entities = ['estn_to_uniprot', 'uniprot_to_seq', 'uniprot_to_domain_seq']

class UniProtAPI:
    """Setting up base_url and headers"""
    def __init__(self, lock):
        self.base_url = "https://rest.uniprot.org"
        self.headers = {"accept": "application/json"}
        self.dicts =  {}
        self.lock = lock

        for entity in entities:
            cache_manager.create_if_not_exists(self.lock, entity)
            self.dicts[entity] = cache_manager.read_whole_json(self.lock, entity)

    def _ensure_cache_entity(self, entity):
        cache_manager.create_if_not_exists(self.lock, entity)
        if entity not in self.dicts:
            self.dicts[entity] = cache_manager.read_whole_json(self.lock, entity)

    """Fetch and return the domain sequence for the given UniProt accession number."""
    def get_domain_seq(self, uniprot_id):
        self._ensure_cache_entity('uniprot_to_domain_seq')
        warnings.warn("Bypassing uniprot: Using curated data from Kannan Lab(TM)")
        return get_domain_seq_from_id(uniprot_id)

        # To do : Check for conflict between curated data and uniprot

        if uniprot_id in self.dicts['uniprot_to_domain_seq'].keys():
            return self.dicts['uniprot_to_domain_seq'][uniprot_id]

        """Fetch protein data from the UniProt API."""
        def get_protein_data(base_url, headers, params):
            response = requests.get(base_url, headers=headers, params=params)
            if not response.ok:
                response.raise_for_status()
            return response.json()

        def extract_sequence_part(data, feature_type='Domain'):
            """Extract the part of the sequence between the start and end positions of a given feature type."""
            sequence = data.get('sequence', {}).get('value', '')
            features = data.get('features', [])
            if(len(features) > 1):
                raise AssertionError("Multiple domains in single sequence")

            # Iterate through features and extract domain sequence
            for feature in features:
                if feature.get('type') == feature_type:
                    start_position = feature.get('location', {}).get('start', {}).get('value', 0)
                    end_position = feature.get('location', {}).get('end', {}).get('value', 0)

                    # Make sure the sequence is long enough to cover the required range
                    if 0 <= start_position < len(sequence) and 0 < end_position <= len(sequence):
                        # Extract the part of the sequence from start to end
                        sequence_part = sequence[start_position - 1:end_position]  # Adjust for 1-based indexing
                        return sequence_part
            return None

        base_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
        params = {
            "fields": [
                "sequence",
                "ft_domain"
            ]
        }
        headers = {
            "accept": "application/json"
        }

        # Get protein data from UniProt API
        try:
            data = get_protein_data(base_url, headers, params)
        except Exception as e:
            return f"Error retrieving data: {e}"

        # Extract the domain sequence
        domain_sequence = extract_sequence_part(data, feature_type='Domain')

        if domain_sequence:
            cache_manager.write_json(self.lock, 'uniprot_to_domain_seq', uniprot_id, domain_sequence)
            return domain_sequence
        else:
            print(uniprot_id)
            raise AssertionError("No domain sequence found")
            return f"No domain feature found for {uniprot_id}."


    """Fetch and return the domain sequence for the given UniProt accession number."""
    def get_seq_from_uniprot_id(self, uniprot_id):
        self._ensure_cache_entity('uniprot_to_seq')
        if uniprot_id in self.dicts['uniprot_to_seq'].keys():
            return self.dicts['uniprot_to_seq'][uniprot_id]

        seq = cache_manager.read_json(self.lock, 'uniprot_to_seq', uniprot_id)
        if seq != None:
            #print("CACHE HIT!")
            self.dicts['uniprot_to_seq'][uniprot_id] = seq
            return seq

        print("CASHE MISS!")
        url = f"{self.base_url}/uniprotkb/{uniprot_id}"
        params = {"fields": ["sequence"]}


        try:
            response = requests.get(url, headers=self.headers, params=params)
            if not response.ok:
                response.raise_for_status()
            data = response.json()
        except Exception as e:
            return f"Error retrieving data: {e}"

        seq = data.get('sequence', {}).get('value', '')
        cache_manager.write_json(self.lock, 'uniprot_to_seq', uniprot_id, seq)

        return seq


    """Orchestrate the job posting and polling to retrieve the mapping."""
    def get_mapping(self, uni_id):

        """Post the mapping job and return the job ID."""
        def post_mapping_job(self, uni_id):
            url = f"{self.base_url}/idmapping/run"
            data = {
                'from': 'UniProtKB_AC-ID',
                'to': 'Gene_Name',
                'ids': uni_id,
                'size' : 1000
            }
            response = requests.post(url, data=data)

            if response.status_code == 200:
                result = response.json()
                job_id = result.get("jobId")
                if job_id:
                    return job_id
                else:
                    return None
            else:
                print(f"Failed to start job: {response.status_code}")
                raise AssertionError(response.text)
                return None

        """Poll the results and return the output once the job is finished."""
        def poll_results(self, job_id):
            status_url = f"{self.base_url}/idmapping/details/{job_id}"
            #print(status_url)
            while True:
                response = requests.get(status_url)

                if response.status_code == 200:
                    status_result = response.json()
                    redirect_url = status_result.get('redirectURL')

                    if redirect_url:
                        result_response = requests.get(redirect_url)

                        if result_response.status_code == 200:
                            final_results = result_response.json()
                            return final_results  # Return the mapping results
                        else:
                            print(f"Failed to fetch final results from redirect URL: {result_response.status_code}")
                            raise AssertionError(result_response.text)
                            break
                    else:
                        print("Mapping still in progress...")
                        time.sleep(30)  # Wait before polling again
                else:
                    print(f"Failed to get job status: {response.status_code}")
                    raise AssertionError(response.text)  # Print any error message
                    break

            return None

        job_id = post_mapping_job(self, uni_id)
        if job_id:
            time.sleep(15)
            return poll_results(self, job_id)
        else:
            raise AssertionError("Mapping from Uniprot_id to Gene_Name not found")

    """Fetch the Gene_name from Uniprot ID"""
    def get_gene_names_from_uniprot_id(self, uni_id):
        # id_mapping has a result limit of 25. Need to chunk it
        batched_uni_ids = [uni_id[i:i + 25] for i in range(0, len(uni_id), 25)]
        gene_names = dict()

        for batch in batched_uni_ids:
            mapping_results = self.get_mapping(batch)

            if "failedIds" in mapping_results:
                for entry in mapping_results["failedIds"]:
                    gene_names[entry] = "Invalid"

            if "results" in mapping_results:
                for entry in mapping_results["results"]:
                    gene_names[entry["from"]] = entry["to"]

        return gene_names

    """Fetch and return the domain sequence for the given UniProt accession number."""
    def get_seq_from_uniprot_id(self, uniprot_id):
        self._ensure_cache_entity('uniprot_to_seq')
        if uniprot_id in self.dicts['uniprot_to_seq'].keys():
            return self.dicts['uniprot_to_seq'][uniprot_id]

        seq = cache_manager.read_json(self.lock, 'uniprot_to_seq', uniprot_id)
        if seq != None:
            #print("CACHE HIT!")
            self.dicts['uniprot_to_seq'][uniprot_id] = seq
            return seq

        print("CASHE MISS!")
        url = f"{self.base_url}/uniprotkb/{uniprot_id}"
        params = {"fields": ["sequence"]}


        try:
            response = requests.get(url, headers=self.headers, params=params)
            if not response.ok:
                response.raise_for_status()
            data = response.json()
        except Exception as e:
            return f"Error retrieving data: {e}"

        seq = data.get('sequence', {}).get('value', '')
        cache_manager.write_json(self.lock, 'uniprot_to_seq', uniprot_id, seq)

        return seq


    """Orchestrate the job posting and polling to retrieve the mapping."""
    def get_uniprot_mapping(self, enst_id):

        """Post the mapping job and return the job ID."""
        def post_mapping_job(self, enst_id):
            url = f"{self.base_url}/idmapping/run"
            data = {
                'from': 'Ensembl_Transcript',
                'to': 'UniProtKB',
                'ids': enst_id,
                'size' : 1000
            }
            response = requests.post(url, data=data)

            if response.status_code == 200:
                result = response.json()
                job_id = result.get("jobId")
                if job_id:
                    return job_id
                else:
                    return None
            else:
                print(f"Failed to start job: {response.status_code}")
                raise AssertionError(response.text)
                return None

        """Poll the results and return the output once the job is finished."""
        def poll_results(self, job_id):
            status_url = f"{self.base_url}/idmapping/details/{job_id}"
            #print(status_url)
            while True:
                response = requests.get(status_url)

                if response.status_code == 200:
                    status_result = response.json()
                    redirect_url = status_result.get('redirectURL')

                    if redirect_url:
                        result_response = requests.get(redirect_url)

                        if result_response.status_code == 200:
                            final_results = result_response.json()
                            return final_results  # Return the mapping results
                        else:
                            print(f"Failed to fetch final results from redirect URL: {result_response.status_code}")
                            raise AssertionError(result_response.text)
                            break
                    else:
                        print("Mapping still in progress...")
                        time.sleep(1000)  # Wait before polling again
                else:
                    print(f"Failed to get job status: {response.status_code}")
                    raise AssertionError(response.text)  # Print any error message
                    break

            return None

        job_id = post_mapping_job(self, enst_id)
        if job_id:
            time.sleep(1)
            return poll_results(self, job_id)
        else:
            raise AssertionError("Mapping from ESTN_id to Uniprot_id not found")

    """Fetch the UniProt ID and validate the sequence."""
    def get_uniprot_ids_from_estn(self, enst_id):
        self._ensure_cache_entity('estn_to_uniprot')
        if len(enst_id) == 1:
            if enst_id[0] in self.dicts['estn_to_uniprot'].keys():
                return {enst_id[0] : self.dicts['estn_to_uniprot'][enst_id[0]]}

            _id = cache_manager.read_json(self.lock, 'estn_to_uniprot', enst_id[0])
            if _id != None :
                self.dicts['estn_to_uniprot'][enst_id[0]] = _id
                return {enst_id[0] : _id}

        # id_mapping has a result limit of 25. Need to chunk it
        batched_enst_ids = [enst_id[i:i + 25] for i in range(0, len(enst_id), 25)]

        uniprot_ids = dict()
        for elem in enst_id:
            if elem in self.dicts['estn_to_uniprot'].keys():
                uniprot_ids[elem] = self.dicts['estn_to_uniprot'][elem]


        for batch in batched_enst_ids:
            #print(batch)
            if all(elem in uniprot_ids for elem in batch):
                #print("Cache hit!!!")
                continue

            print("CASHE MISS!")
            mapping_results = self.get_uniprot_mapping(batch)
            #print(mapping_results)

            if "failedIds" in mapping_results:
                for entry in mapping_results["failedIds"]:
                    uniprot_ids[entry] = "Invalid"
                    cache_manager.write_json(self.lock, 'estn_to_uniprot', entry, "Invalid")


            if "results" in mapping_results:
                for entry in mapping_results["results"]:
                    uniprot_id = entry["to"]["primaryAccession"]
                    uniprot_ids[entry["from"]] = uniprot_id
                    seq_from_mapping = entry["to"]["sequence"]["value"]
                    #seq_from_id = self.get_seq_from_uniprot_id(uniprot_id)

                    cache_manager.write_json(self.lock, 'estn_to_uniprot', entry["from"], uniprot_id)
        return uniprot_ids

# Example usage
if __name__ == "__main__":

    api = UniProtAPI(multiprocessing.Lock())
    print(api.get_gene_names_from_uniprot_id(["A0A2I3JL00"])["A0A2I3JL00"])

    #sys.exit(0)

    #Q8IXL6
    #print(api.get_domain_seq("Q13627"))
    #print(api.get_domain_seq("Q13627"))
    #sys.exit(0)

    #for entry in get_well_studied_human_kinase_ids():
    #    print(api.get_domain_seq(entry))
    #sys.exit(0)


    #for entity in entities:
    #    print(f"Length of {entity} dictionary = {len(api.dicts[entity])}")
    #sys.exit(0)

    uniprot_id = "Q9Y261"
    print(api.get_seq_from_uniprot_id(uniprot_id))
    #sys.exit(0)

    #enst_id = ["ENST00000221978", "ENST00000346753", "ENST00000371263", "ENST00000372108", "ENST00000470487"]
    enst_id = ["ENST00000311015"]
    try:
        uniprot_id = api.get_uniprot_ids_from_estn(enst_id)
        print(f"UniProt ID: {uniprot_id}")
    except Exception as e:
        print(f"Error: {e}")
