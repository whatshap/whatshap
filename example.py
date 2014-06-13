#!/usr/bin/env python3
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
		#print(repr(line))
		s = ''
		for pos, c in enumerate(line, 1):
			if c == ' ':
				continue
			print('{} {} {} {} : '.format(pos*10, bits[c], c, quality), end='', file=f)
		print('# 55 : NA', file=f)

output = subprocess.check_output(['build/dp', '--all_het', 'example.wif'], shell=False).decode()
print(output)

result = output.split('\n')[0]

bases = ''
for field in result.split(' : ')[:-2]:
	pos, base, bit, qual = field.split(' ')
	assert base == bit
	print(pos, base, qual)
	bases += base
print(bases)
#72 T 0 39 : 84 T 0 41 : # 26 : NA
