from __future__ import print_function
import sys

from eip1829.eip1829 import eip1829, EIP1829Error


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
			is_error = line[0]
			expected = tuple(line[1:3])
			arguments = line[3:]
			try:
				actual = eip1829(*arguments)
			except EIP1829Error as ex:
				if not is_error:
					raise ex
			else:
				if is_error:
					print("Arguments:", arguments)
					print("Expected:", expected)
					print("Actual:", actual)
					raise RuntimeError("Expected error, got no error!")
			if expected != actual:
				print("Arguments:", arguments)
				print("Expected:", expected)
				print("Actual:", actual)
				raise RuntimeError("Test case mismatch")

	return 0


if __name__ == "__main__":
	sys.exit(main(sys.argv))
