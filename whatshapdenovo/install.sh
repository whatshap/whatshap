
basedir=`pwd`
rm -rf bin
cd $basedir/tools/HaploConduct
make
cd $basedir/tools/CONSENT
./install.sh
# cd $basedir/bin
# ln -fs ../tools/CONSENT/CONSENT-correct .
cd $basedir

