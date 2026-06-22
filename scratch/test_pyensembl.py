import pyensembl
ensembl = pyensembl.EnsemblRelease(release=111, species="human")
# Download/install database if not already done
ensembl.download()
ensembl.index()

db = ensembl.db
print("Database type:", type(db))
print("Database path:", db.database_path)

# Let's try to query the exon table
try:
    conn = db.connect()
    cursor = conn.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cursor.fetchall()
    print("Tables in Ensembl DB:", tables)
    
    # Check structure of 'exon' table
    cursor.execute("PRAGMA table_info(exon);")
    exon_info = cursor.fetchall()
    print("Exon table columns:", exon_info)
    
    # Check if we can select gene_id, start, end from exon
    cursor.execute("SELECT gene_id, start, end FROM exon LIMIT 5;")
    print("Exon sample rows:", cursor.fetchall())
except Exception as e:
    print("Error querying database:", e)
