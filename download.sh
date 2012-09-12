#!/bin/bash
# need root 
# get depend data and program from internet
# and update the original ChiLin config
# at the same time check all dependency
help()
{
    cat <<HELP
This is the help for install ChiLin and all dependency file
USAGE EXAMPLE: sh install.sh -p pathd
pathd: where to store download data
HELP
exit 0
}

while [ -n "$1" ]; do
    case $1 in
		-h)help;shift 1;;
		-p)pathd=$2; shift 2;;
		--)shift; break;; # options end
		-*) echo "error: no such option $1. -h for help";exit 1;;
		*) break;;
		esac
done

# choose where to install
PS3="Choose (1-2)"
echo "Choose from the list below."
select name in global personal
do
    break
done

if [ $name = global ]; then
    echo "install to global"
else
    echo "get ChiLin data to $pathd"
fi

if [ -d $pathd ]
then
    echo "directory exists"
else
    mkdir $pathd
fi
cd $pathd

# for judging machine type for download the right
# source or binary
mach() {
    MACH="set |grep -a i386"
    mac64="set | (grep x86_64 && grep apple)"
    mac32="set | (grep i386 && grep apple)"
    linux= "set | (grep i386 && grep linux)"
    if [ -z echo $mac64 ]; then
			echo "Your machine type is 64 mac machine"
		elif [ -z echo $mac32 ]; then
			echo "Your machine type is 32 mac machine"
    elif [ -z echo $linux ]; then
			echo "Your machine type is linux machine"
    fi
}

get() {
	file=$1
	cmd="type -a $2"
	if ! $cmd > /dev/null
		then
		    if ! wget --version >/dev/null 2>/dev/null ; then
					if ! curl --version >/dev/null 2>/dev/null ; then
								    echo "Please install wget or curl somewhere in your PATH"
											    exit 1
							fi
							curl -o `basename $1` -L $1
									return $?
			    else
					wget $1 -O `basename $1`
							return $?
			    fi
		    else
		    echo "$2 already installed"
	fi
}

install() {
    suffix=$1
    tar xvfz $2   # for samtools and bowtie
	  cd $2
    if [ ${suffix##*.} = py ];then
	    python setup.py install
	    return $?
    elif [ ${suffix#*.} = cc ]; then
	    make
	    cp $2 /usr/bin/
	    return $?
	  else
			cp $3 /usr/bin/
			return $?
    fi
    cd -
}

# fastqc
if get http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.10.1.zip fastqc;
then
    install java fastqc_* fastqc
fi

# bowtie
species="\
mm9 \
hg19 \
"
BowtieIndexBase="ftp://ftp.cbcb.umd.edu/pub/data/bowtie_indexes/"
for index in $species; do
    if [ ! -f ${index}.ebwt.zip ]; then
			if get ${BowtieIndexBase}/${index}.ebwt.zip bowtie; then
				    echo "Bowtie index download successfully"
			fi
		fi
done

CURL_BASE="\
http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.8/bowtie-0.12.8-src.zip/download \
http://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2/download \
http://bedtools.googlecode.com/files/BEDTools.v2.16.2.tar.gz \
"
for url in $CURL_BASE:
do
    if get "http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.8/bowtie-0.12.8-src.zip/download" bowtie
    then
        install cc bowtie-* bowtie
    fi
done

if get http://cistrome.org/~hanfei/conservation/hg19/vertebrate/ conservation_plot.py
then
    echo "conservation preparation ready"
fi

binary_base="\
http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedClip \
http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig \
"

if type -a macs2
then
   echo "macs2 already install"
else:
   wget --no-check-certificate https://github.com/taoliu/MACS/zipball/master -O macs2.tar.gz
   echo "macs2 downloaded"
   tar xvfz macs2.tar.gz
   cd taoliu-MACS*
   python setup.py install
fi

if type -a ceas
then
   echo "ceas already installed"
else
   wget -c http://liulab.dfci.harvard.edu/CEAS/src/CEAS-Package-1.0.2.tar.gz
   tar xvfz CEAS-Package-1.0.2.tar.gz
   cd CEAS*
   # refgene
   wget -c http://liulab.dfci.harvard.edu/CEAS/src/mm9.refGene.gz
   wget -c http://liulab.dfci.harvard.edu/CEAS/src/hg19.refGene.gz
   gunzip *.gz
fi

motif
if type -a MDSeqPos.py
then
   echo "MDSeqPos already installed"
else
   wget -c https://bitbucket.org/cistrome/cistrome-applications-harvard/get/e82ed15a486b.tar.gz -O cistrome
   tar xvfz *.tar.gz
fi

ChiLin source
get https://bitbucket.org/shenglinmei/chilin/get/6e472d7e69a6.tar.gz ChiLin py
tar xvfz *.tar.gz
cd shenglinmei*
python setup.py install

if $?
then
   echo "ChiLin data config succeed!:)"
fi
