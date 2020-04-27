# Notes

## Required Binaries For Utility Scripts

#### get_dbSNPvcf_contig_mappings.sh
The phoenix_resources.ini currently points to binaries available in /home/jkeats today for NCBI eUtils and JSON.awk

```
### NCBI eUtils
## Install instructions: https://www.ncbi.nlm.nih.gov/books/NBK179288/
cd ~
/bin/bash
perl -MNet::FTP -e \
'$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1);
$ftp->login; $ftp->binary;
$ftp->get("/entrez/entrezdirect/edirect.tar.gz");'
gunzip -c edirect.tar.gz | tar xf -
rm edirect.tar.gz
builtin exit
export PATH=${PATH}:$HOME/edirect >& /dev/null || setenv PATH "${PATH}:$HOME/edirect"
./edirect/setup.sh

### JSON.awk
cd ~
wget https://github.com/step-/JSON.awk/archive/1.3.tar.gz
tar xvzf 1.3.tar.gz
```