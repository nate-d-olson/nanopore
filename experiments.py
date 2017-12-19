import numpy
import pandas
import collections

METADATA_URL = "https://docs.google.com/spreadsheets/d/1dmhJKHVHEnZdELOTm9jl0uevGlArZSs9zR97oRY_IcY/gviz/tq?tqx=out:csv"

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
                "completed": istrue(row.Completed)
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