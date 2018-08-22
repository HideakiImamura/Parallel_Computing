import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


data = pd.read_csv("results.csv")
data_mat = np.zeros((8, 9), dtype=np.float)
for _, x in data.iterrows():
    data_mat[int(x['nthread'] - 1), int(np.log10(int(x['size'])))-1] = float(x['par_time'])
# print(data_mat)

########################################
for i, y in enumerate(data_mat.T):
    x = np.arange(0, len(y)) + 1
    plt.plot(x, y, label="$i$ = " + str(i))
plt.legend()
plt.xlabel("# of thread")
plt.ylabel("Elapsed time [mm]")
plt.title("Scalability: size of array is $10^i$")
plt.savefig("scalability.pdf")
plt.clf()

for i, y in enumerate(data_mat.T):
    x = np.arange(0, len(y)) + 1
    plt.plot(x, np.log(y), label="$i$ = " + str(i))
plt.legend()
plt.xlabel("# of thread")
plt.ylabel("$\log$(Elapsed time) [mm]")
plt.title("Scalability: size of array is $10^i$")
plt.savefig("scalability_log.pdf")
plt.clf()

########################################

ps = []
for i in range(0, 9):
    s = 8
    a = data_mat[0][i] / data_mat[7][i]
    p = (1 / a - 1) / (1 / s - 1)
    print(p)
    ps.append(p)
plt.plot(np.arange(5)+4, ps[4:])
plt.title("Proportion of parallelizable")
plt.xlabel("$i$: size of array is $10^i$")
plt.ylabel("Proportion of parallelizable")
plt.savefig("parallel.pdf")
