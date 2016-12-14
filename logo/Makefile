all: whatshap_logo.pdf whatshap_logo.png whatshap_logo_small.png

%.pdf: %.svg
	inkscape -A $@ $<

%.png: %.svg
	inkscape -d 160 -e $@ $<

