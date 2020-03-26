from matplotlib import pyplot as plt
from random_var import RandomVariable

rv = RandomVariable([0.7, 0.29, 0.01])
seqs = rv.get_all_sequences(2)
wts, comp = rv.partition_weakly_typical_set(seqs, 0.9)
print('wts:', wts)  # {(0, 1), (1, 0), (0, 0), (1, 1)}
print('comp: wts', comp)  # {(1, 2), (2, 1), (2, 0), (2, 2), (0, 2)}

sts, comp = rv.partition_strongly_typical_set(seqs, 0.9)
print('sts:', sts)
print('comp sts:', comp, '\n')


# The weakly typical set contains most of the probability mass even
# though it does not usually contain the most probable sequence
rv = RandomVariable([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.4])
seqs = rv.get_all_sequences(1)  # {(0,), (1,), (2,), (3,), (4,), (5,)}
wts, comp = rv.partition_weakly_typical_set(seqs, 0.9)
print('wts:', wts)
print('comp wts:', comp)

# According to the Weak Asymptotic Equipartition Property, the proportion of sequences
# in the weakly epsilon-typical set should be > 1-epsilon as n → infinity.
rv = RandomVariable([0.6, 0.4])
ratios = [(lambda a, b: len(a)/(len(a)+len(b)+0.01))(*rv.partition_weakly_typical_set(rv.get_all_sequences(i), 0.2)) for i in range(1, 15)]
# ['0.0', '0.4987', '0.7490', '0.6246', '0.7810', '0.8748', '0.9296', '0.8515', '0.9081', '0.9443', '0.9667', '0.9267', '0.9537', '0.9712']
plt.plot(ratios)
plt.show()

# According to the Strong Asymptotic Equipartition Property, the proportion of sequences
# in the strongly delta-typical set should be > 1-delta as n → infinity.
rv = RandomVariable([0.6, 0.4])
ratios = [(lambda a, b: len(a)/(len(a)+len(b)+0.01))(*rv.partition_strongly_typical_set(rv.get_all_sequences(i), 0.2)) for i in range(1, 15)]
plt.plot(ratios)
plt.show()