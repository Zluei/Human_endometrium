import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("str_all_count.csv");
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("str_all_count.loom",x.X.transpose(),row_attrs,col_attrs);
