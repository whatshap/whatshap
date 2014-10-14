#!/usr/bin/env python3
"""
Example for paired-end reads.

The correct output would be two connected components.
"""

import subprocess

IN = """
0011         1101
 1100        00101
  1001        01010
        101"""
#1234567890123456789
#111001  101  001010
s = IN.strip()

bits = { '0': 'A', '1': 'C', 'E': 'G' }

with open('example.wif', 'w') as f:
	quality = 22
	for line in s.split('\n'):
		s = ''
		for pos, c in enumerate(line, 1):
			if c == ' ':
				continue
			print('{} {} {} {} : '.format(pos*10, bits[c], c, quality), end='', file=f)
		print('# 55 : NA', file=f)

print('Creating an "example.wif" file with the following content:')
with open('example.wif') as f:
	print(f.read())

print('Running "dp" program.')
output = subprocess.check_output(['build/dp', '--all_het', 'example.wif'], shell=False).decode()
print('Output by "dp" program:')
print(output, end='')
print('End of output.\n')

print('Result (position, base, quality):')
result = output.split('\n')[0]

bases = ''
for field in result.split(' : ')[:-2]:
	pos, base, bit, qual = field.split(' ')
	assert base == bit
	print(pos, base, qual)
	bases += base

print('\nResult as bitstring:', bases)
#72 T 0 39 : 84 T 0 41 : # 26 : NA
