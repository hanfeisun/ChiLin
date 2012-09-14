#!/bin/bash
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
              need root gcc numpy ghostscript and django
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

cp chilin/lib/db/*.bed $bin
cp chilin/lib/db/*.txt $bin
unzip chilin/lib/db/DHS.zip -d $bin


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
        machine="linux.x86_64"
	echo "Your machine type is linux machine"
    else
        machine="linux.i386"
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
install other bedClip
install other bedGraphToBigWig
unzip fastqc_v0.10.1.zip -d $bin
cd $bin
chmod 755 fastqc
cd -

## for sourceforge special
wget -c http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.8/bowtie-0.12.8-src.zip/download -O bowtie.zip
wget -c http://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2/download -O samtools.tar.bz2
unzip bowtie.zip
install cc bowtie.zip
tar xvfj samtools.tar.bz2
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
					unzip ${BowtieIndexBase}/${index}.ebwt.zip -d $bin
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
#   ceas, conservation_plot, mdseqpos, bigwiggle_correlation.py
#   venn_diagram.py
# included data
#   refgene
if [ ! -f cistrome-app.tar.gz ]; then
    wget --no-check-certificate \
        https://bitbucket.org/cistrome/cistrome-applications-harvard/get/e82ed15a486b.tar.gz -O cistrome-app.tar.gz
    tar xvfz cistrome-app.tar.gz
    cd cistrome-applications-harvard*
    # ceas
    cd published-packages/CEAS
    python setup.py install --prefix=`dirname $bin`
    cd gdb && gunzip hg19* && gunzip mm9*
    cp hg19* mm9* $bin
    cd ../../cistrome-extra-apps
    python setup.py install --prefix=`dirname $bin`
    # mdseqpos problem
    export CC=gcc  # upper case
    cd ../../mdseqpos
    python setup.py install --prefix=`dirname $bin`
fi
#
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
		# vertebrate
		wget -c http://cistrome.org/~hanfei/conservation/$s/vertebrate/${c}.bw -O conservation_$s/${c}.bw
    done
    cp -r conservation_$s $bin
done

cd - #back to the chilin directory

echo "\
[UserInfo]
User = 
datasetID = 
species = {{ species }}
factor = 
treatpath = 
controlpath = 
OutputDirectory = 
# should add fastq or bam suffix
[bowtie]
BAMTOFASTQ = $bin/bam2fastx
BOWTIE_MAIN = $bin/bowtie
BOWTIE_GENOME_INDEX_PATH = $bin/BowtieIndex/{{ species }}
nBOWTIE_MAX_ALIGNMENT = 1

[samtools]
SAMTOOLS_MAIN = $bin/samtools
SAMTOOLS_CHROM_LEN_PATH = $bin/chromInfo_{{ species }}.txt
CHROM_LEN_BED_PATH = $bin/chr_limit_{{ species }}.bed

[macs]
MACS_MAIN = $bin/macs2
BEDGRAPHTOBIGWIG_MAIN = $bin/bedGraphToBigWig
bedclip = $bin/bedClip

[bedtools]
INTERSECTBED_MAIN = $bin/intersectBed

[bed2bam]
BEDTOBAM_MAIN = $bin/bedToBam

[ceas]
CEAS_MAIN = $bin/ceasBW
ceas_ex = $bin/ceas-exon
CEAS_GENETABLE_PATH = $bin/{{ species }}.refGene
CHROM_LEN = $bin/chromInfo_{{ species }}.txt
CEAS_PROMOTER_SIZES = 
CEAS_BIPROMOTER_SIZES = 
CEAS_REL_DIST = 

[conservation]
CONSERV_PLOT_MAIN = $bin/conservation_plot.py
PEAKS_NUM = 3000
CONSERV_PLOT_PHAST_PATH = $bin/conservation_{{ species }}
# use vertebrate Phascon

[correlation]
WIG_CORRELATION_MAIN = $bin/wig_correlation_bigwig_input_only.py
WIG_CORRELATION_STEP = 10
WIG_CORRELATION_METHOD = mean 
WIG_CORRELATION_MIN = 2
WIG_CORRELATION_MAX = 50

[venn]
VENN_DIAGRAM_MAIN = $bin/venn_diagram.py
DHS_BED_PATH = $bin/DHS/DHS_{{ species }}.bed 
VELCRO_PATH = $bin/wgEncodeHg19ConsensusSignalArtifactRegions.bed
# velcro only available for hg19

[seqpos]
SEQPOS_MAIN = $bin/MDSeqPos.py
SEQPOS_TOP_PEAKS = 1000
SEQPOS_MDSCAN_WIDTH = 200
SEQPOS_MDSCAN_TOP_PEAKS = 200
SEQPOS_MDSCAN_TOP_PEAKS_REFINE = 500
SEQPOS_WIDTH = 600
SEQPOS_PVALUE_CUTOFF = 0.001
SEQPOS_MOTIF_DB_SELECTION = cistrome.xml

[QC]
FASTQC_MAIN = $bin/fastqc
FILTERDUP_SPECIES = {{ filterdup_species }}
" > ChiLin.conf.new
if $?
then
   echo "ChiLin data config succeed!:)"
fi
