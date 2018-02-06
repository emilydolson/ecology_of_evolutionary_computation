import glob
import sys
import pandas as pd


frames = []

for res_dir in glob.glob(sys.argv[1]):
    fit_data = pd.read_csv(res_dir+"/fitness.csv", index_col="update")
    sys_data = pd.read_csv(res_dir+"/systematics.csv", index_col="update")
    ex_sys_data = pd.read_csv(res_dir+"/extra_systematics.csv", index_col="update")
    sel_data = pd.read_csv(res_dir+"/selection_info.csv", index_col="update")
    pop_data = pd.read_csv(res_dir+"/population.csv", index_col="update")

    all_data = pd.concat([fit_data, sys_data, ex_sys_data, sel_data, pop_data],axis=1)
    all_data["seed"] = res_dir.split("_")[-1]
    all_data["selection"] = res_dir.split("_")[0]
    if len(res_dir.split("_")) == 3:
        all_data["sel_param"] = res_dir.split("_")[1]
    else:
        all_data["sel_param"] = -1

    frames.append(all_data)

all_data = pd.concat(frames)

all_data.to_csv("all_data.csv")
