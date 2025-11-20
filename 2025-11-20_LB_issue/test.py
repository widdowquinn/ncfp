from Bio import Entrez

Entrez.email = "ANother@example.com"

dbname = "protein"
linkdbname = "protein_nuccore"
query_id = "EPD1771162.1"
query_id = "XP_004520832.1"

# Current ncfp implementation
elink = Entrez.elink(dbfrom=dbname, linkname=linkdbname, id=query_id, retmode="text")
print(f"\n{elink=}")

print(f"\n{dir(elink)=}")

try:
    matches = Entrez.read(elink)
    print(f"\n{matches=}")
except RuntimeError as err:
    print(f"\n{err=}")

# Using "ipg" as at https://www.biostars.org/p/234693/
ipg_data = Entrez.efetch(db=dbname, rettype="ipg", retmode="xml", id=query_id)
print(f"\n{ipg_data=}")

ipg_matches = Entrez.read(ipg_data)
print(f"\n{ipg_matches=}")
