import numpy
import pandas
import collections

METADATA_URL = "https://docs.google.com/spreadsheets/d/1dmhJKHVHEnZdELOTm9jl0uevGlArZSs9zR97oRY_IcY/export?gid=0&format=csv"
METADATA_URL_NOTT = "https://docs.google.com/spreadsheets/d/1dmhJKHVHEnZdELOTm9jl0uevGlArZSs9zR97oRY_IcY/export?gid=908644433&format=csv"

def istrue(val):
    try:
        if numpy.isnan(val): return False
    except TypeError:
        pass
    if str(val).lower() in ["y", "yes", "true", "1"]:
        return True
    if str(val).lower() in ["n", "no", "false", "0"]:
        return False
    raise Exception(f"is this true? {val}")
    

def load_experiment_metadata():
    experiments = _load_our_experiments()
    experiments.update(_load_their_experiments())
    return experiments
    
def _load_their_experiments():
    df = pandas.read_csv(METADATA_URL_NOTT, parse_dates=["Date"])
    print(df)
    indices = collections.defaultdict(dict)
    for flowcell_id, group in df.groupby("Flowcell ID"):
        assert len(group["Type"].unique()) == 1
        
        indices[flowcell_id]["flowcell_type"] = list(group["Type"])[0]
        indices[flowcell_id]["datasets"] = []
        
        for _, row in group.iterrows():
            if pandas.isnull(row.Run):
                continue
            cur_run = {
                "date":row.Date,
                "name": row.Run,
                "sample": row["Sample ID"],
                "basecalled": True,
            }
            indices[flowcell_id]["datasets"].append(cur_run)
            
        if len(indices[flowcell_id]["datasets"]) == 0:
            del indices[flowcell_id]

    # sort by date
    indices = collections.OrderedDict((k,indices[k]) for k in 
                                      sorted(indices, reverse=True,
                                             key=lambda x: max(d["date"] for d in indices[x]["datasets"])))
    # pprint.pprint(indices)

    return indices

    
def _load_our_experiments():
    df = pandas.read_csv(METADATA_URL, parse_dates=["Date"])
        
    indices = collections.defaultdict(dict)
    for flowcell_id, group in df.groupby("Flowcell ID"):
        assert len(group["Type"].unique()) == 1
        
        indices[flowcell_id]["flowcell_type"] = list(group["Type"])[0]
        indices[flowcell_id]["datasets"] = []
        
        for _, row in group.iterrows():
            if pandas.isnull(row.Run):
                continue
            cur_run = {
                "date":row.Date,
                "kit": row.Kit,
                "name": row.Run,
                "completed": istrue(row.Completed),
                "sample": row["Sample ID"],
                "basecalled": False,
            }
            indices[flowcell_id]["datasets"].append(cur_run)
            
        if len(indices[flowcell_id]["datasets"]) == 0:
            del indices[flowcell_id]

    # sort by date
    indices = collections.OrderedDict((k,indices[k]) for k in 
                                      sorted(indices, reverse=True,
                                             key=lambda x: max(d["date"] for d in indices[x]["datasets"])))
    # pprint.pprint(indices)

    return indices


if __name__ == '__main__':
    import pprint
    pprint.pprint(load_experiment_metadata())
