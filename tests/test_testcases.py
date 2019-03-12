import sys

from eip1829.eip1829 import eip1829


def main(argv):
	filename = 'tests/testcases.csv'
	if len(argv) >= 2:
		filename = argv[1]
	with open(filename, 'r') as handle:
		first = True
		for line in handle:
			if first:
				first = False
				continue
			line = [int(_, 16) for _ in line.strip().split(',')]
			if len(line) < 8:
				raise RuntimeError('Invalid number of arguments')
			expected = tuple(line[:2])
			arguments = line[2:]
			actual = eip1829(*arguments)
			if expected != actual:
				print("Arguments:", arguments)
				print("Expected:", expected)
				print("Actual:", actual)
				raise RuntimeError("Test case mismatch")

	return 0


if __name__ == "__main__":
	sys.exit(main(sys.argv))
