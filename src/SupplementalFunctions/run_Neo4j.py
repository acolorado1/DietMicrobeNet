# script developed for Neo4j integration
# better for graphs >1M nodes, so far they have been smaller 
# can replace run_graph.py if needed 
# requires neo4j desktop 

from neo4j import GraphDatabase
import pandas as pd
from tqdm import tqdm
import argparse as arg

# === Define Neo4j write functions ===
def clear_database(tx):
    """Delete all nodes and relationships."""
    tx.run("MATCH (n) DETACH DELETE n")
def create_nodes(tx, df):
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Creating nodes", ncols=80):
        tx.run("""
            MERGE (p:Compound {c_id: $compound})
            SET p.origin = $origin,
                p.assoc_food = $assoc_food,
                p.freq = $freq
        """, compound=row["compound"],
             origin=row["origin"],
             assoc_food=row["assoc_food"],
             freq=row["freq"])

def create_edges(tx, df):
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Creating edges", ncols=80):
        tx.run("""
            MATCH (c1:Compound {c_id: $compound1})
            MATCH (c2:Compound {c_id: $compound2})
            MERGE (c1)-[r:BECOMES]->(c2)
            SET 
                r.reaction = $reaction,
                r.KOs = $KOs,
                r.organisms = $organisms, 
                r.abundance = $abundance
        """, compound1=row["compound1"],
             compound2=row["compound2"],
             reaction=row["reaction"],
             KOs=row["KOs"],
             organisms=row["organisms"],
             abundance=row["abundance"])

# === Pattern queries ===
cypher1 = """
MATCH (c1:Compound {origin:"food"})-[r:BECOMES]->(c2:Compound {origin:"microbe"})
RETURN 
    c1.c_id AS compound1_id,
    c1.origin AS compound1_origin,
    c1.assoc_food AS compound1_assoc_food,
    c1.freq AS compound1_freq,
    r.reaction AS reaction,
    r.KOs AS KOs,
    r.organisms AS organisms,
    r.abundance AS abundance,
    c2.c_id AS compound2_id,
    c2.origin AS compound2_origin,
    c2.assoc_food AS compound2_assoc_food,
    c2.freq AS compound2_freq
"""
cypher2 = """
MATCH (c1:Compound {origin:"food"})-[r:BECOMES]->(c2:Compound {origin:"both"})
RETURN 
    c1.c_id AS compound1_id,
    c1.origin AS compound1_origin,
    c1.assoc_food AS compound1_assoc_food,
    c1.freq AS compound1_freq,
    r.reaction AS reaction,
    r.KOs AS KOs,
    r.organisms AS organisms,
    r.abundance AS abundance,
    c2.c_id AS compound2_id,
    c2.origin AS compound2_origin,
    c2.assoc_food AS compound2_assoc_food,
    c2.freq AS compound2_freq
"""

cypher3 = """
MATCH (c1:Compound {origin:"both"})-[r:BECOMES]->(c2:Compound {origin:"both"})
RETURN 
    c1.c_id AS compound1_id,
    c1.origin AS compound1_origin,
    c1.assoc_food AS compound1_assoc_food,
    c1.freq AS compound1_freq,
    r.reaction AS reaction,
    r.KOs AS KOs,
    r.organisms AS organisms,
    r.abundance AS abundance,
    c2.c_id AS compound2_id,
    c2.origin AS compound2_origin,
    c2.assoc_food AS compound2_assoc_food,
    c2.freq AS compound2_freq
"""

# === define main ===
def main(): 
    parser = arg.ArgumentParser(description='Connect to Neo4j, create and query graph')
    parser.add_argument('--n', type=str, required=True, help='node file path')
    parser.add_argument('--e', type=str, required=True, help='edge file path')
    parser.add_argument('--uri', type=str, required=True, help='Neo4j URI instance')
    parser.add_argument('--user', type=str, required=True, help='Neo4j username for instance')
    parser.add_argument('--p', type=str, required=True, help='Neo4j password for instance')
    parser.add_argument('--o', type=str, required=True, help='output file path')
    args = parser.parse_args()

    # === Create Neo4j driver ===
    driver = GraphDatabase.driver(
        args.uri,
        auth=(args.user, args.p)
    )
    driver.verify_connectivity()
    print("✅ Connected to Neo4j.\n")

    # === Load CSVs ===
    nodes_df = pd.read_csv(args.n)
    edges_df = pd.read_csv(args.e)
    print(f"Loaded {len(nodes_df)} nodes and {len(edges_df)} edges from CSV files.\n")
    
    # === Run everything ===
    with driver.session() as session:
        print("🧹 Clearing existing graph data...")
        session.execute_write(clear_database)

        print("🧩 Creating graph in Neo4j...\n")
        session.execute_write(create_nodes, nodes_df)
        session.execute_write(create_edges, edges_df)

        print("\n🔍 Running pattern query 1: food -> microbe...")
        result1 = session.run(cypher1)
        df1 = pd.DataFrame([r.data() for r in result1])

        print("\n🔍 Running pattern query 2: food -> both...")
        result2 = session.run(cypher2)
        df2 = pd.DataFrame([r.data() for r in result2])

        print("\n🔍 Running pattern query 3: both -> both...")
        result3 = session.run(cypher3)
        df3 = pd.DataFrame([r.data() for r in result3])

    driver.close()

    # === Output results ===
    print(f"\n✅ Queries complete. Retrieved {len(df1)+len(df2)} relationships.")
    df = pd.concat([df1, df2, df3], 
                   ignore_index=True,
                   sort=False)

    # === Save results ===
    df.to_csv(args.o, index=False)
    print(f"\n📁 Saved results to: {args.o}")

if __name__ == "__main__":
    main()