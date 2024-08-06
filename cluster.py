from pandas_ods_reader import read_ods
import pandas as pd


def parse_address(address):
    """Convert hierarchical address to a list of integers."""
    return list(map(int, address.split('.')))

def group_snp_addresses(snp_list, threshold):
    """Group SNP addresses by their hierarchical similarity up to a given threshold level."""
    # Parse the SNP addresses into a list of lists
    parsed_addresses = [parse_address(snp) for snp in snp_list]
    
    # Dictionary to map hierarchical prefix to group number
    prefix_to_group = {}
    group_id = 0
    
    # List to store the group number for each SNP address
    groups = []
    
    for address in parsed_addresses:
        # Build the prefix from the address up to the threshold level
        prefix = tuple(address[:threshold])
        if prefix in prefix_to_group:
            assigned_group = prefix_to_group[prefix]
        else:
            # New prefix, assign a new group number
            assigned_group = group_id
            prefix_to_group[prefix] = group_id
            group_id += 1
        
        groups.append(assigned_group)
    
    return groups



def cluster(start_date = '2014-01-01', end_date = '2015-12-31'):
    start_date = pd.to_datetime(start_date, format='%Y-%m-%d')
    end_date = pd.to_datetime(end_date, format='%Y-%m-%d')
    path = "metadata.ods"
    df = read_ods(path)
    df['RECEIPT_datetime'] = pd.to_datetime(df['RECEIPT'], format='%Y-%m-%d')
    print(df)
    df['Cluster'] = group_snp_addresses(df["SNP"].values.tolist(), 3)
    df = df.sort_values("RECEIPT_datetime")
    df = df.sort_values("Cluster")
    print(df)
    df.to_csv("meta_with_cluster.csv", index=False)

    df_subset = df[(df['RECEIPT_datetime'] >= start_date) & (df['RECEIPT_datetime'] <= end_date)]
    isolates = df_subset["SRA Accession"].values.tolist()

    dfs = [group for _, group in df_subset.groupby(["Cluster"])]
    cluster_list = []
    for df_ in dfs:
        df_cluster = pd.DataFrame(df_["SRA Accession"].values.tolist())
        culster_id = df_["Cluster"].values.tolist()[0]
        filepath = f"cluster_{culster_id}.csv"
        print(filepath)
        cluster_list.append(filepath)
        df_cluster.to_csv(filepath, index=False)
    
    return isolates, cluster_list
