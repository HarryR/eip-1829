m = 21888242871839275222246405745257275088696311157297823662689037894645226208583
alpha = 0
beta = 3
g1 = (1, 2)

g1_times_two = (1368015179489954701390400359078579693043519447331113978918064868415326638035, 9918110051302171585080402603319702774565515993150576347155970296011118125764)
g1_times_three = (3353031288059533942658390886683067124040920775575537747144343083137631628272, 19321533766552368860946552437480515441416830039777911637913418824951667761761)
g1_times_five = (10744596414106452074759370245733544594153395043370666422502510773307029471145, 848677436511517736191562425154572367705380862894644942948681172815252343932)

from eip1829.eip1829 import eip1829

assert eip1829(m, alpha, beta, [2], [g1]) == g1_times_two
assert eip1829(m, alpha, beta, [3], [g1]) == g1_times_three
assert eip1829(m, alpha, beta, [5], [g1]) == g1_times_five