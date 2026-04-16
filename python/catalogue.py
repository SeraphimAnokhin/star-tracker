
import numpy as np

data = np.genfromtxt('hipparcos-voidmain.csv', delimiter=',', skip_header=1).transpose()

m = data[5]
RA = data[8] * np.pi / 180
dec = data[9] * np.pi / 180


RA = RA[m < 5]
dec = dec[m < 5]
m = m[m < 5]

dec = np.delete(dec, np.isnan(RA))
m = np.delete(m, np.isnan(RA))
RA = np.delete(RA, np.isnan(RA))

n = (m.size)

dist_mat = np.zeros(n * (n - 1) // 2)
pair_table = np.zeros(n * (n - 1) // 2, dtype=[('k', np.int32), ('l', np.int32), ('dist', np.float64)])

for i in range(n):
    for j in range(i + 1, n):
        #print(i, j, np.sin(dec[i]) * np.sin(dec[j]) + np.cos(dec[i]) * np.cos(dec[j]) * np.cos(RA[i] - RA[j]))
        dist_mat[(int)(i * (n - 1.5 - 0.5 * i) + j - 1)] = np.arccos(np.sin(dec[i]) * np.sin(dec[j]) + np.cos(dec[i]) * np.cos(dec[j]) * np.cos(RA[i] - RA[j])) * 180 / np.pi
        pair_table[(int)(i * (n - 1.5 - 0.5 * i) + j - 1)][0] = i
        pair_table[(int)(i * (n - 1.5 - 0.5 * i) + j - 1)][1] = j
        pair_table[(int)(i * (n - 1.5 - 0.5 * i) + j - 1)][2] = dist_mat[(int)(i * (n - 1.5 - 0.5 * i) + j - 1)]

#print(pair_table.shape)
pair_table.sort(0, order='dist')
#print(pair_table.shape)

delta_hash = 0.05
max_dist = 60
hash_table = np.zeros((int)(max_dist / delta_hash), dtype=np.int32)

hash_table[0] = 0
for i in range(1, (int)(max_dist / delta_hash)):
    hash_table[i] = hash_table[i - 1]
    #print(i, pair_table[hash_table[i]][2])
    while (pair_table[hash_table[i]][2] < i * delta_hash):
        hash_table[i] += 1
        #print(i, pair_table[hash_table[i]][2])

np.savetxt("hash_table.csv", hash_table, fmt="%d")
np.savetxt("pair_table.csv", pair_table[pair_table['dist'] <= max_dist], fmt="%d;%d;%.5lf")
np.savetxt("dist_mat.csv", dist_mat, fmt="%.5lf")
np.savetxt("star_coord.csv", np.array((RA, dec)).T * 180 / np.pi, fmt="%.5lf;%.5lf")
np.savetxt("star_mag.csv", m, fmt="%.5lf")
