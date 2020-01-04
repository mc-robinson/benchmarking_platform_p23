# Working with the Data Folder

This folder usually contains 4 large files:

- `chembl_smiles_df.pkl`, 131M
- `dud_smiles_df.pkl`, 536M
- `dud_missing_smiles_df.pkl`, 35M
- `muv_smiles_df.pkl`, 736M

These four files provide all of the necessary CAS fingerprints for the benchmark. However, they all are too big to store on GitHub. The files can be downloaded once inside the `data` directory from a remote storage server as follows:

`wget https://cas-benchmark-data.s3.eu-west-2.amazonaws.com/CAS-github-data-files/chembl_smiles_df.pkl`
`wget https://cas-benchmark-data.s3.eu-west-2.amazonaws.com/CAS-github-data-files/dud_missing_smiles_df.pkl`
`wget https://cas-benchmark-data.s3.eu-west-2.amazonaws.com/CAS-github-data-files/dud_smiles_df.pkl`
`wget https://cas-benchmark-data.s3.eu-west-2.amazonaws.com/CAS-github-data-files/muv_smiles_df.pkl`

After downloading, they should all be in the `data` directory before using the platform.






