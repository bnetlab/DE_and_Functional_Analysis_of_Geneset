# geoQuery

library(GEOquery)
gse <- getGEO("GSE21653", GSEMatrix = TRUE)
show(gse)

gpl97 <- getGEO('GPL97')
Meta(gpl97)$title
head(Meta(gpl97)$series_id)

# download fastq from ebi

system("wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/ERR458791/ERR458791.fastq.gz")
system("wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/ERR459179/ERR459179.fastq.gz")
ERR459179

require(XML)
data <- xmlParse("ena.xml")
xml_data <- xmlToList(data)

for (i in 1:672){
print(xml_data[i]$RUN$IDENTIFIERS$PRIMARY_ID)
}

system("wget http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=ERR458793&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes")