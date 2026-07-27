import urllib.request
import json
import re

url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=GSE120575[Accession]&retmode=json"
with urllib.request.urlopen(url) as response:
    data = json.loads(response.read().decode())
    idlist = data.get("esearchresult", {}).get("idlist", [])
    if idlist:
        print("GDS IDs:", idlist)

import ftplib
try:
    ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd('/geo/series/GSE120nnn/GSE120575/suppl')
    files = ftp.nlst()
    print("Files in FTP:")
    for f in files:
        print(f)
    ftp.quit()
except Exception as e:
    print("FTP error:", e)
