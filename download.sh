#!/bin/bash
# get depend data and program from internet
# and update the original ChiLin config
# at the same time check all dependency
#------------------------
help()
{
    cat <<HELP
usage: $0 install [-p DIR]
optional arguments:
  -p DIR		specify the directory to store the downloaded data DEFAULT: DOWNLOADED
HELP
}

# options
while [ -n "$1" ]; do
    case $1 in
        -h)help;shift 1;;
        -p)pathd=$2; shift 2;;
        --)shift; break;; # options end
        -*) echo "error: no such option $1. -h for help";exit 1;;
        *) break;; 
    esac
done

if [ -n $pathd ]
then
    echo "input path right"
else
    echo "no such directory"
    help
    exit 0
fi
mkdir -p $pathd

## choose for all users or single user
PS3="Choose (1-2)"
echo "Choose from the list below."
select name in global personal
do
   break
done

if [ $name = global ]; then
    u=`whoami`
fi

if [[ $u = root ]]
then
    echo "install for all users"
else
    echo "No root authority, please ask your administrator!"
    echo "Or try personal install, use local directory for bin and data"
fi

# users could customize executive bin path
read -p "choose where the bin is " bin
read -p "choose where the data is" data
if [ ! -n $bin ];then
    if (echo $PATH|grep $bin) ; then echo 1;
    else
        echo "input bin not in $PATH, use default /usr/bin"
        bin=/usr/bin
    fi
fi
if [ ! -f $data ]; then
    mkdir $data
fi

# cp built-in data
cp chilin/db/*.bed $data
cp chilin/db/*.txt $data
unzip chilin/db/DHS.zip -d $data

# change into pathd to download and install
cd $pathd

# for judging machine type for download the right
# source or binary
if uname -a | grep Darwin; then
    machine="macOSX.i386"
    echo "your machine is mac"
elif uname -a |grep Linux; then
    if uname -a | grep x86; then
        machine="linux.x86_64"
        echo "Your machine type is linux.x86_64"
    else
        machine="Your machine type is linux.i386"
    fi
fi

get() {
    file=$1
    if ! wget --version >/dev/null 2>/dev/null ; then
        if ! curl --version >/dev/null 2>/dev/null ; then
            echo "Please install wget or curl somewhere in your PATH"
            exit 1
        fi
        curl -o $2 -L $1
        return $?
    else
        wget $1 -O $2
        return $?
    fi
}

install() {
    suffix=$1
    if [ ${suffix##*.} = py ];then
        if python -V 2>&1 | grep 2.7; then
            echo "installining $2"
            python setup.py install --prefix=`dirname $bin`
        else
            echo "Please install python 2.7 version"
            return $?
        fi
    elif [ ${suffix#*.} = cc ]; then
        name=`echo $2|cut -d . -f 1`
        cd $name-*
        make
        if [ -d bin ]; then
            cp -r bin/* $bin
        else
            cp $name $bin
        fi
        cd -
        return $?
    else
        cp `basename $2` $bin
        return $?
    fi
}

# depend data
# bowtie index
# TODO
# for ab solid colorspace alignment use color index hgc
species="\
mm9 \
hg19 \
"
BowtieIndexBase="ftp://ftp.cbcb.umd.edu/pub/data/bowtie_indexes/"
for index in $species; do
    if [ ! -f ${index}.ebwt.zip ]; then
        if get ${BowtieIndexBase}/${index}.ebwt.zip ${index}.ebwt.zip; then
            echo "Bowtie index download successfully"
            unzip ${BowtieIndexBase}/${index}.ebwt.zip -d $data
        fi
    fi
done
# depend software
URL_BASE="\
http://bedtools.googlecode.com/files/BEDTools.v2.16.2.tar.gz \
http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.10.1.zip \
http://hgdownload.cse.ucsc.edu/admin/exe/$machine/bedClip \
http://hgdownload.cse.ucsc.edu/admin/exe/$machine/bedGraphToBigWig \
"
for i in $URL_BASE;do
    if [ ! -f `basename $i` ]; then
        if get $i `basename $i`; then
            echo "download $i ready"
        fi
    fi
done
# sourceforge special
if [ ! -f bowtie.zip ]; then
    get http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.8/bowtie-0.12.8-src.zip/download bowtie.zip
    if [ ! -f samtools.tar.bz2 ]; then
        get http://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2/download  samtools.tar.bz2
    fi
fi

# install
install other bedClip
install other bedGraphToBigWig
unzip fastqc_v0.10.1.zip -d $bin >/dev/null
cd $bin
chmod 755 $bin/FastQC/fastqc
cd -

unzip bowtie.zip >/dev/null
echo "install bowtie"
install cc bowtie.zip
tar xvfj samtools.tar.bz2 >/dev/null 2>&1
install cc samtools.tar.bz2
echo "install samtools"

# github certificate macs2
if type -a macs2
then
    echo "macs2 already install"
else
    wget --no-check-certificate https://github.com/taoliu/MACS/zipball/master -O taoliu-MACS2.tar.gz
    echo "Latest macs2 downloaded from github"
    export CC=gcc  # upper case
    tar xvfz taoliu-MACS2.tar.gz>/dev/null
    cd taoliu*
    install py macs2
fi

# mercurial cistrome-applications
# included program
#   ceas, conservation_plot, mdseqpos, bigwiggle_correlation.py
#   venn_diagram.py
# included data
#   refgene
if [ ! -f cistrome-app.tar.gz ]; then
    wget --no-check-certificate \
        https://bitbucket.org/cistrome/cistrome-applications-harvard/get/e82ed15a486b.tar.gz -O cistrome-app.tar.gz
fi
tar xvfz cistrome-app.tar.gz >/dev/null 2>&1
cd cistrome-cistrome-applications-harvard*

# ceas
cd published-packages/CEAS
install py CEAS
cd gdb && gunzip *.refGene.gz
cp -rf *.refGene $bin
cd ../../../cistrome-extra-apps built-ins
install py built-ins

# mdseqpos 
export CC=gcc  # upper case
cd ../mdseqpos
cp lib/settings.py.example lib/settings.py
install py mdseqpos

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
for s in $species;do
    for c in $BASE_CHRS;do
        mkdir conservation_$s 2>/dev/null
        wget -c http://cistrome.org/~hanfei/conservation/$s/vertebrate/${c}.bw -O conservation_$s/${c}.bw
    done
    cp -r conservation_$s $data
done

echo "\
[basis]
user = foo
id = 2012
time = 2102xxx
species = {{ species }}
factor = ESR1
treat = absolutepath
control = absolutepath
output = absolutepath

[bowtie]
path = $bin/bowtie
# for color sequence using color index
index = $data/{{ species }}
bam2fq = $bin/bamToFastq
maxalign = 1

[samtools]
path = $bin/samtools
chrlen = $data/chromInfo_{{ species }}.txt
chrbed = $data/chr_limit_{{ species }}.bed

[macs]
path = $bin/macs2
bg2bw = $bin/bedGraphToBigWig
bedclip = $bin/bedClip

[bedtools]
path = $bin/intersectBed
dhs = $data/DHS_{{ species }}.bed
velcro = $data/wgEncodeHg19ConsensusSignalArtifactRegions.bed

[ceas]
path = $bin/ceasBW
exon = $bin/ceas-exon
refgene = $data/{{ species }}.refGene
chrlen = $data/chromInfo_{{ species }}.txt
promoter_sizes =
bipromoter_sizes =
rel_dist =

[conservation]
path = $bin/conservation_plot.py
peaks = 3000
width = 4000
phast = $data/conservation_{{ species }}

[correlation]
path = $bin/bigwig_correlation.py
wig_correlation_step = 10
wig_correlation_method = mean
wig_correlation_min = 2
wig_correlation_max = 50

[venn]
path = $bin/venn_diagram.py

[seqpos]
path = $bin/MDSeqPos.py
species = {{ filterdup_species }}
summitsnumber = 1000
mdscan_width = 200
mdscan_top_peaks = 200
seqpos_mdscan_top_peaks_refine = 500
seqpos_width = 600
seqpos_pvalue_cutoff = 0.001
db = cistrome.xml

[qc]
path = $bin/FastQC/fastqc
species = {{ filterdup_species }}
" > chilin/template/ChiLin.conf.sample
if [ $? -eq 0 ]; then
    echo "Now automatic  generated conf is chilin/template/ChiLinjinja.conf.sample"
    echo "template conf has been replaced by the new conf"
    cp chilin/template/ChiLin.conf.sample chilin/conf/ChiLin.conf
fi

# install Chilin
cd ../../.. #back to the chilin root directory
echo "========================================="
if [ $? -eq 0 ]; then
    echo "ChiLin data config succeed!:)"
fi
echo "Please run ChiLin setup.py"


