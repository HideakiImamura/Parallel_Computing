import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


data = pd.read_csv("result.csv")
data_mat = np.zeros((8, 10), dtype=np.float)
for _, x in data.iterrows():
    data_mat[int(x['nthread'] - 1), int(np.log2(int(x['height']) // 60))] = float(x['time'])
# print(data_mat)

########################################
for i, y in enumerate(data_mat.T):
    x = np.arange(0, len(y)) + 1
    plt.plot(x, y, label="$i$ = " + str(i))
plt.legend()
plt.xlabel("# of thread")
plt.ylabel("Elapsed time [mm]")
plt.title("Scalability: size of image is $2^i$" + r"$\times(60,100)$")
plt.savefig("scalability.pdf")
plt.clf()

for i, y in enumerate(data_mat.T):
    x = np.arange(0, len(y)) + 1
    plt.plot(x, np.log(y), label="$i$ = " + str(i))
plt.legend()
plt.xlabel("# of thread")
plt.ylabel("$\log$(Elapsed time) [mm]")
plt.title("Scalability: size of image is $2^i$" + r"$\times(60,100)$")
plt.savefig("scalability_log.pdf")
plt.clf()

########################################

ps = []
for i in range(0, 10):
    s = 8
    a = data_mat[0][i] / data_mat[7][i]
    p = (1 / a - 1) / (1 / s - 1)
    print(p)
    ps.append(p)
plt.plot(np.arange(10), ps)
plt.title("Proportion of parallelizable")
plt.xlabel("$i$: size of image is $2^i$" + r"$\times(60,100)$")
plt.ylabel("Proportion of parallelizable")
plt.savefig("parallel.pdf")
