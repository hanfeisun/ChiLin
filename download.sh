#!/bin/bash
# need root 
# get depend data and program from internet
# and update the original ChiLin config
# at the same time check all dependency

#------------------------
# # chilin source
# get https://bitbucket.org/shenglinmei/chilin/get/6e472d7e69a6.tar.gz
# tar xvfz *.tar.gz
# cd shenglinmei*
# python setup.py install



help()
{
    cat <<HELP
Instruction: download all the depend and install automatically
USAGE EXAMPLE: bash install.sh -p pathd
    -p pathd: where to store download data
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

if [ -z $pathd ]; then
    help
    echo "please input the path for containing download data"
fi

# choose for all users or single user
PS3="Choose (1-2)"
echo "Choose from the list below."
select name in global personal
do
    break
done
if [ $name = global ]; then
    u=`whoami`
    if [ $u = root ]; then
        echo "install for all users"
    else
        echo "need root authority, ask your administrator!"
        exit 1
    fi
else
    echo "get ChiLin data to $pathd for single user"
fi

# Users could customize executive bin path
read -p "choose where the bin is " bin
if [ ! -n $bin ];then
    echo "use default bin path /usr/bin"
    bin=/usr/bin
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
if uname -a | grep Darwin; then
    if uname -a | grep x86; then
        machine="macOSX.ppc"
	echo "Your machine type is 64 mac machine"
    else
        machine="macOSX.i386"
        echo "Your machine type is 32 mac machine"
    fi
elif uname -a |grep Linux; then
    if uname -a | grep x86; then
        machine="Linux.x86_64"
	echo "Your machine type is linux machine"
    else
        machine="Linux.i386"
    fi
fi

get() {
    file=$1
    # cmd="type -a $2"
    # if ! $cmd > /dev/null
    # then
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
    # else
    #     echo "$2 already installed"
    # fi
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
        name=`echo $2|cut -d . -f 1`
	cp $name $bin
	return $?
    else
	cp `basename $2` $bin
	return $?
    fi
    cd -
}

# depend source code and binary
# machine type specific 
URL_BASE="\
http://bedtools.googlecode.com/files/BEDTools.v2.16.2.tar.gz \
http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.10.1.zip \
http://hgdownload.cse.ucsc.edu/admin/exe/$machine/bedClip \
http://hgdownload.cse.ucsc.edu/admin/exe/$machine/bedGraphToBigWig \
"
for i in $URL_BASE;do
    if get $i; then
        echo "download $i ready"
    fi
done

# for sourceforge special
wget -c http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.8/bowtie-0.12.8-src.zip/download -O bowtie.zip
wget -c http://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2/download -O samtools.tar.bz2
install cc bowtie.zip
install cc samtools.tar.bz2

# bowtie index
# you may add more index
species="\
mm9 \
hg19 \
"
BowtieIndexBase="ftp://ftp.cbcb.umd.edu/pub/data/bowtie_indexes/"
for index in $species; do
    if [ ! -f ${index}.ebwt.zip ]; then
	if get ${BowtieIndexBase}/${index}.ebwt.zip; then
	    echo "Bowtie index download successfully"
	fi
    fi
done

# github certificate macs2
if type -a macs2
then
   echo "macs2 already install"
else:
   wget --no-check-certificate https://github.com/taoliu/MACS/zipball/master -O macs2.tar.gz
   echo "Latest macs2 downloaded from github"
   install py macs2.tar.gz
fi
# mercurial cistrome-applications
# included program
#   ceas, conservation_plot
# included data
#   refgene
if [ ! -f cistrome-app.tar.gz ]; then
    wget --no-check-certificate \
        https://bitbucket.org/cistrome/cistrome-applications-harvard/get/e82ed15a486b.tar.gz -O cistrome-app.tar.gz
    tar xvfz cistrome-app.tar.gz
    cd cistrome-applications-harvard
    # ceas
    cd published-packages/CEAS
    python setup.py install
    cd gdb && gunzip hg19* && gunzip mm9*
    cp hg19* mm9* $bin
    cd ../../cistrome-extra-apps
    python setup.py install
    # mdseqpos problem
fi
# special download need copy sub directory
# conservation Phascon bigwiggle
BASE_CHRS="\
chr1 \
chr2 \
chr3 \
chr4 \
chr5 \
chr6 \
chr7 \
chr8 \
chr9 \
chr10 \
chr11 \
chr12 \
chr13 \
chr14 \
chr15 \
chr16 \
chr17 \
chr18 \
chr19 \
chr20 \
chr21 \
chr22 \
chrX \
chrY \
chrM"
for c in $BASE_CHRS;do
    for s in $species;do
        mkdir conservationhg19
        wget -c http://cistrome.org/~hanfei/conservation/$s/vertebrate/${c}.bw -O conservationhg19/${c}.bw
    done
done

#wget -c -r -np -nd  http://cistrome.org/~hanfei/conservation/hg19/vertebrate/ # will get index.html

# update config for default
# ChiLin.conf.new
echo 


if $?
then
   echo "ChiLin data config succeed!:)"
fi
