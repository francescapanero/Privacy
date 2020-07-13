import numpy as np
import pandas as pd

n = 1000

zipf_param = (1.0526,1.1765, 1.3333, 1.5385,1.8182,2.2222,2.8571,4,6.6667,20)

frequencies = np.unique(np.random.zipf(zipf_param[0],n),return_counts=True)

import pandas as pd 
pd.DataFrame(np_array).to_csv("path/to/file.csv")
