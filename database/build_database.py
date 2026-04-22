import sqlite3
import pandas as pd

db_path = "database/neoepitope.db"
csv_path = "data/final_tables/final_prioritized_peptides.csv"

conn = sqlite3.connect(db_path)

df = pd.read_csv(csv_path)
df.to_sql("final_candidates", conn, if_exists="replace", index=False)

conn.close()
print(f"Database created: {db_path}")
