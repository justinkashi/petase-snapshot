INSTALL INTERPRO (more details at https://interproscan-docs.readthedocs.io/en/v5/HowToDownload.html)

mkdir ~/interproscan
cd ~/interproscan

# download
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.76-107.0/interproscan-5.76-107.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.76-107.0/interproscan-5.76-107.0-64-bit.tar.gz.md5

# verify checksum
md5sum -c interproscan-5.76-107.0-64-bit.tar.gz.md5

# extract
tar -xzf interproscan-5.76-107.0-64-bit.tar.gz
cd interproscan-5.76-107.0

# setup
python3 setup.py -f interproscan.properties
chmod +x interproscan.sh# INSTALL INTERPRO (visit https://interproscan-docs.readthedocs.io/en/v5/HowToDownload.html for more details) 


