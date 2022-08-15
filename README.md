# **<span style="color:#4086f7;font-size:20px">Gene and Primer Sequence Analysis for SARS-CoV-2, EGFR(Non Small Lung Cancer Cell), Influenza DNAs</span>**

<a href="https://www.canarydetect.com/"><img src="https://raw.githubusercontent.com/Owaiskhan9654/Gene-Sequence-Primer-/main/Canary%20Logo.jpg"></a>

<center><b>Visit Company website by clicking on the Logo</b></center>


# **<span style="color:#4086f7;font-size:20px">Gene and Primer Sequence Analysis for SARS-CoV-2, EGFR(Non Small Lung Cancer Cell), Influenza DNAs</span>**

 ### This Notebook is Created by [**Owais Ahmad**](https://www.linkedin.com/in/owaiskhan9654/) for SARS-CoV-2 Gene Sequence Analysis.
 ### This is for Medical Research Purpose only.
 
 
 
- **Email owaiskhan9654@gmail.com**
- **Contact +919515884381**

## First time setup - Generating credentials
1. Sign in or create a new IDT account [here](https://www.idtdna.com/site/account?returnurl=/site/account/api).
2. Go to your user name drop down menu at the top right of the page, and select My account.
3. Click the API access link.
4. Click the Request new API key button.
5. Append those 4 secret in Kaggle Add-ons
6. Comment down for any help.

### How can I check my Oligo primers to ensure there are no significant primer design issues?
- The difference between melting temperatures (Tm) of the primers should be less than 5°C.
- The GC content should be between 35-80% or equivalent to the product being amplified.
- The Delta G value of any self-dimers, hairpins, and heterodimers should be weaker (more positive) than -9.0 kcal/mole. Positive numbers indicate that the actual secondary structure shown will not form at all.
- Avoid 3' complementarity between the two primers to prevent primer dimers.
The IDT OligoAnalyzer APIs can be used to assess these different criteria for a proposed oligo. 
#### [Reference](https://sg.idtdna.com/pages/support/faqs/how-can-i-check-my-pcr-primers-using-the-oligoanalyzer-program-to-ensure-there-are-no-significant-primer-design-issues-)

![dna-spiral-genetics-twisted-wallpaper-preview.jpg](attachment:ee1e1a3c-d9f3-41a0-8183-7e07ebe5eba9.jpg)

# **<span style="color:#4086f7;">Basic Library Installations</span>**


```python
!pip install -q openpyxl
```

    [33mWARNING: Running pip as the 'root' user can result in broken permissions and conflicting behaviour with the system package manager. It is recommended to use a virtual environment instead: https://pip.pypa.io/warnings/venv[0m[33m

    [0m


```python
from __future__ import print_function
from base64 import b64encode
import json
from urllib import request, parse
import pandas as pd
import requests
import pandas as pd
from tqdm.notebook import tqdm
from functools import reduce
from kaggle_secrets import UserSecretsClient
```

# **<span style="color:#4086f7; font-size:18px">Secret Token Fetch from IDT-DNA. Note this will expire in 10 Mins and if you have a very long set of sequence to perform analysis then call below function again to get new Secret token</span>**


```python
def get_bearer_token(client_id, client_secret, idt_username, idt_password):


    authorization_string = b64encode(bytes(client_id + ":" + client_secret, "utf-8")).decode()
    request_headers = { "Content-Type" : "application/x-www-form-urlencoded",
                        "Authorization" : "Basic " + authorization_string }
                    
    data_dict = {   "grant_type" : "password", "scope" : "test","username" : idt_username,"password" : idt_password }
    request_data = parse.urlencode(data_dict).encode()

    post_request = request.Request("https://www.idtdna.com/Identityserver/connect/token", data = request_data, headers = request_headers,method = "POST")

    response = request.urlopen(post_request)

    body = response.read().decode()
    
    if (response.status != 200):
        raise RuntimeError("Request failed with error code:" + str(response.status) + "\nBody:\n" + body)
    
    body_dict = json.loads(body)
    return body_dict["access_token"]
    
user_secrets = UserSecretsClient()
secret_value_0 = user_secrets.get_secret("IDT_password_here")
secret_value_1 = user_secrets.get_secret("IDT_username_here")
secret_value_2 = user_secrets.get_secret("client_id_here")
secret_value_3 = user_secrets.get_secret("client_secret_here")
    
client_id = str(secret_value_2)
client_secret = str(secret_value_3)
idt_username = str(secret_value_1)
idt_password = str(secret_value_0)


token = get_bearer_token(client_id, client_secret, idt_username, idt_password)
print("Secret Token Fetched (This will expire in 10 minutes): ",token)

```

    Secret Token Fetched (This will expire in 10 minutes):  c7adab9d36e0f8a43664dc398cc78322
    

# **<span style="color:#4086f7; font-size:18px">Here SARS-CoV-2 Gene Sequences are fetched from csv file which is stored on GitHub. This we get from blasting on NCBI</span>**


```python
url="https://github.com/Owaiskhan9654/Gene-Sequence-Primer-/blob/main/NEB%20Primer%20Sequence.xlsx?raw=true"
# Sequence fetch for Analysis

response = requests.get(url)

dest = 'GENE Primer Sequence.xlsx'

with open(dest, 'wb') as file:
    file.write(response.content)

GENE_df = pd.read_excel("GENE Primer Sequence.xlsx", sheet_name=2,header=1).dropna()
GENE_df.reset_index(drop=True,inplace=True)
GENE_df.columns
```




    Index(['Primer name', 'Sequence', 'Synthesis scale'], dtype='object')




```python
GENE_df.to_csv('SARS-CoV-2_Primer_Sequences.csv',index=False)
```


```python
GENE_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Primer name</th>
      <th>Sequence</th>
      <th>Synthesis scale</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Gene-E1-F3</td>
      <td>TGAGTACGAACTTATGTACTCAT</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Gene-E1-B3</td>
      <td>TTCAGATTTTTAACACGAGAGT</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Gene-E1-FIP</td>
      <td>ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Gene-E1-BIP</td>
      <td>TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Gene-E1-LF</td>
      <td>CGCTATTAACTATTAACG</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>5</th>
      <td>Gene-E1-LB</td>
      <td>GCGCTTCGATTGTGTGCGT</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>6</th>
      <td>N2-F3</td>
      <td>ACCAGGAACTAATCAGACAAG</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>7</th>
      <td>N2-B3</td>
      <td>GACTTGATCTTTGAAATTTGGATCT</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>8</th>
      <td>N2-FIP</td>
      <td>TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>9</th>
      <td>N2-BIP</td>
      <td>CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>10</th>
      <td>N2-LF</td>
      <td>GGGGGCAAATTGTGCAATTTG</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>11</th>
      <td>N2-LB</td>
      <td>CTTCGGGAACGTGGTTGACC</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>12</th>
      <td>ACTB-F3</td>
      <td>AGTACCCCATCGAGCACG</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>13</th>
      <td>ACTB-B3</td>
      <td>AGCCTGGATAGCAACGTACA</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>14</th>
      <td>ACTB-FIP</td>
      <td>GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>15</th>
      <td>ACTB-BIP</td>
      <td>CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>16</th>
      <td>ACTB-LF</td>
      <td>TGTGGTGCCAGATTTTCTCCA</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>17</th>
      <td>ACTB-LB</td>
      <td>CGAGAAGATGACCCAGATCATGT</td>
      <td>10 nm</td>
    </tr>
  </tbody>
</table>
</div>



![image.png](attachment:6c728c92-b0bf-4b85-b01e-a8383e7945e6.png)


```python
Sequence_dict = {}

for i in GENE_df.index:
    Sequence_dict[GENE_df['Sequence'][i]] = GENE_df['Primer name'][i]

Sequence_dict

Primer_dict = {}

for i in GENE_df.index:
    Primer_dict[GENE_df['Primer name'][i]] = GENE_df['Sequence'][i]

Primer_dict
```




    {'Gene-E1-F3': 'TGAGTACGAACTTATGTACTCAT',
     'Gene-E1-B3': 'TTCAGATTTTTAACACGAGAGT',
     'Gene-E1-FIP': 'ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG',
     'Gene-E1-BIP': 'TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT',
     'Gene-E1-LF': 'CGCTATTAACTATTAACG',
     'Gene-E1-LB': 'GCGCTTCGATTGTGTGCGT',
     'N2-F3': 'ACCAGGAACTAATCAGACAAG',
     'N2-B3': 'GACTTGATCTTTGAAATTTGGATCT',
     'N2-FIP': 'TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC',
     'N2-BIP': 'CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA',
     'N2-LF': 'GGGGGCAAATTGTGCAATTTG',
     'N2-LB': 'CTTCGGGAACGTGGTTGACC',
     'ACTB-F3': 'AGTACCCCATCGAGCACG',
     'ACTB-B3': 'AGCCTGGATAGCAACGTACA',
     'ACTB-FIP': 'GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA',
     'ACTB-BIP': 'CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC',
     'ACTB-LF': 'TGTGGTGCCAGATTTTCTCCA',
     'ACTB-LB': 'CGAGAAGATGACCCAGATCATGT'}




```python
%%time

headers = {
    'Content-Type': 'application/json',
    'Accept': 'application/json',
    'Authorization': 'Bearer '+token,
}

Primer_Name = []
Sequence = []
Complement = []
length = []
GCContent = []
MeltTemp = []
NmoleOD = []
OligoConc_list = []
MgConc_list = []
NaConc_list = []
dNTPsConc_list = []
NucleotideType_list=[]
count = 0
for i in tqdm(GENE_df.index):
    Primer_Name1=GENE_df.iloc[i]['Primer name']
#     print(Primer_Name1)
    oligo_Conc_dict={'F3':0.2,'B3':0.2,'Fip':1.6,'Bip':1.6,'Lf':0.4,'Bf':0.4} 
    dntp_Conc_dict={'F3':1.4,'B3':1.4,'Fip':1.4,'Bip':1.4,'Lf':1.4,'Bf':1.4} 
    mg_Conc_dict={'F3':8,'B3':8,'Fip':8,'Bip':8,'Lf':8,'Bf':8} 
    Na_Conc_dict={'F3':50,'B3':50,'Fip':50,'Bip':50,'Lf':50,'Bf':50}
    if 'F3' in Primer_Name1:
        oligo_Conc,dntp_Conc,mg_Conc,Na_Conc=oligo_Conc_dict['F3'],dntp_Conc_dict['F3'],mg_Conc_dict['F3'],Na_Conc_dict['F3']
    elif 'B3' in Primer_Name1:
        oligo_Conc,dntp_Conc,mg_Conc,Na_Conc=oligo_Conc_dict['B3'],dntp_Conc_dict['B3'],mg_Conc_dict['B3'],Na_Conc_dict['B3']
    elif 'FIP' in Primer_Name1:
        oligo_Conc,dntp_Conc,mg_Conc,Na_Conc=oligo_Conc_dict['Fip'],dntp_Conc_dict['Fip'],mg_Conc_dict['Fip'],Na_Conc_dict['Fip']
    elif 'BIP' in Primer_Name1:
        oligo_Conc,dntp_Conc,mg_Conc,Na_Conc=oligo_Conc_dict['Bip'],dntp_Conc_dict['Bip'],mg_Conc_dict['Bip'],Na_Conc_dict['Bip']
    elif 'LF' in Primer_Name1:
        oligo_Conc,dntp_Conc,mg_Conc,Na_Conc=oligo_Conc_dict['Lf'],dntp_Conc_dict['Lf'],mg_Conc_dict['Lf'],Na_Conc_dict['Lf']
    elif 'BF' in Primer_Name1:
        oligo_Conc,dntp_Conc,mg_Conc,Na_Conc=oligo_Conc_dict['Bf'],dntp_Conc_dict['Bf'],mg_Conc_dict['Bf'],Na_Conc_dict['Bf']
    
    #print(oligo_Conc)
    NucleotideType="DNA"
    data = '{  "Sequence": "' + GENE_df.iloc[i]['Sequence'] + '",  "NaConc": '+  str(Na_Conc)+\
    ',  "MgConc": '+  str(mg_Conc)+',   "dNTPsConc": '+  str(dntp_Conc)+',  "OligoConc": '+  str(oligo_Conc)+\
    ',   "NucleotideType": "DNA" }'
    response = requests.post(
        'https://www.idtdna.com/Restapi/v1/OligoAnalyzer/Analyze',
        headers=headers,
        data=data)
    json_data = json.loads(response.text)
    Primer_Name.append(Sequence_dict[json_data['Sequence'].replace(" ", '')])
    Sequence.append(json_data['Sequence'])
    Complement.append(json_data['Complement'])
    length.append(json_data['Length'])
    GCContent.append(json_data['GCContent'])
    MeltTemp.append(json_data['MeltTemp'])
    NmoleOD.append(json_data['NmoleOD'])
    OligoConc_list.append(json_data['OligoConc'])
    MgConc_list.append(mg_Conc)
    NaConc_list.append(Na_Conc)
    dNTPsConc_list.append(dntp_Conc)
    NucleotideType_list.append(NucleotideType)


df1=pd.DataFrame({"Primer Name":Primer_Name, "Sequence":Sequence,"Complement":Complement,"OligoConc":OligoConc_list,"Na+ Conc":NaConc_list,\
                  "Mg++ Conc":MgConc_list,"dNTPs Conc":dNTPsConc_list,"Nucleotide Type":NucleotideType_list,"length":length,"GCContent":GCContent,"MeltTemp":MeltTemp,"NmoleOD":NmoleOD,})


```


      0%|          | 0/18 [00:00<?, ?it/s]


    CPU times: user 419 ms, sys: 25.9 ms, total: 445 ms
    Wall time: 18.2 s
    

# **<span style="color:#4086f7; font-size:22px">In this below DataFrame you can check for Length of Formation, Its GC content, Melting Temperature, and NmoleOD  </span>**


```python
df1
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Primer Name</th>
      <th>Sequence</th>
      <th>Complement</th>
      <th>OligoConc</th>
      <th>Na+ Conc</th>
      <th>Mg++ Conc</th>
      <th>dNTPs Conc</th>
      <th>Nucleotide Type</th>
      <th>length</th>
      <th>GCContent</th>
      <th>MeltTemp</th>
      <th>NmoleOD</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Gene-E1-F3</td>
      <td>TGA GTA CGA ACT TAT GTA CTC AT</td>
      <td>ATG AGT ACA TAA GTT CGT ACT CA</td>
      <td>0.2</td>
      <td>50</td>
      <td>8</td>
      <td>1.4</td>
      <td>DNA</td>
      <td>23</td>
      <td>34.8</td>
      <td>61.0</td>
      <td>4.40</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Gene-E1-B3</td>
      <td>TTC AGA TTT TTA ACA CGA GAG T</td>
      <td>ACT CTC GTG TTA AAA ATC TGA A</td>
      <td>0.2</td>
      <td>50</td>
      <td>8</td>
      <td>1.4</td>
      <td>DNA</td>
      <td>22</td>
      <td>31.8</td>
      <td>60.4</td>
      <td>4.57</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Gene-E1-FIP</td>
      <td>ACC ACG AAA GCA AGA AAA AGA AGT TCG TTT CGG AA...</td>
      <td>CTG TCT CTT CCG AAA CGA ACT TCT TTT TCT TGC TT...</td>
      <td>1.6</td>
      <td>50</td>
      <td>8</td>
      <td>1.4</td>
      <td>DNA</td>
      <td>42</td>
      <td>42.9</td>
      <td>75.8</td>
      <td>2.26</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Gene-E1-BIP</td>
      <td>TTG CTA GTT ACA CTA GCC ATC CTT AGG TTT TAC AA...</td>
      <td>ACG TGA GTC TTG TAA AAC CTA AGG ATG GCT AGT GT...</td>
      <td>1.6</td>
      <td>50</td>
      <td>8</td>
      <td>1.4</td>
      <td>DNA</td>
      <td>44</td>
      <td>40.9</td>
      <td>75.4</td>
      <td>2.41</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Gene-E1-LF</td>
      <td>CGC TAT TAA CTA TTA ACG</td>
      <td>CGT TAA TAG TTA ATA GCG</td>
      <td>0.4</td>
      <td>50</td>
      <td>8</td>
      <td>1.4</td>
      <td>DNA</td>
      <td>18</td>
      <td>33.3</td>
      <td>53.4</td>
      <td>5.69</td>
    </tr>
    <tr>
      <th>5</th>
      <td>Gene-E1-LB</td>
      <td>GCG CTT CGA TTG TGT GCG T</td>
      <td>ACG CAC ACA ATC GAA GCG C</td>
      <td>0.4</td>
      <td>50</td>
      <td>8</td>
      <td>1.4</td>
      <td>DNA</td>
      <td>19</td>
      <td>57.9</td>
      <td>68.6</td>
      <td>5.86</td>
    </tr>
    <tr>
      <th>6</th>
      <td>N2-F3</td>
      <td>ACC AGG AAC TAA TCA GAC AAG</td>
      <td>CTT GTC TGA TTA GTT CCT GGT</td>
      <td>0.2</td>
      <td>50</td>
      <td>8</td>
      <td>1.4</td>
      <td>DNA</td>
      <td>21</td>
      <td>42.9</td>
      <td>61.2</td>
      <td>4.52</td>
    </tr>
    <tr>
      <th>7</th>
      <td>N2-B3</td>
      <td>GAC TTG ATC TTT GAA ATT TGG ATC T</td>
      <td>AGA TCC AAA TTT CAA AGA TCA AGT C</td>
      <td>0.2</td>
      <td>50</td>
      <td>8</td>
      <td>1.4</td>
      <td>DNA</td>
      <td>25</td>
      <td>32.0</td>
      <td>62.2</td>
      <td>4.20</td>
    </tr>
    <tr>
      <th>8</th>
      <td>N2-FIP</td>
      <td>TTC CGA AGA ACG CTG AAG CGG AAC TGA TTA CAA AC...</td>
      <td>GGC CAA TGT TTG TAA TCA GTT CCG CTT CAG CGT TC...</td>
      <td>1.6</td>
      <td>50</td>
      <td>8</td>
      <td>1.4</td>
      <td>DNA</td>
      <td>42</td>
      <td>47.6</td>
      <td>77.6</td>
      <td>2.43</td>
    </tr>
    <tr>
      <th>9</th>
      <td>N2-BIP</td>
      <td>CGC ATT GGC ATG GAA GTC ACA ATT TGA TGG CAC CT...</td>
      <td>TAC ACA GGT GCC ATC AAA TTG TGA CTT CCA TGC CA...</td>
      <td>1.6</td>
      <td>50</td>
      <td>8</td>
      <td>1.4</td>
      <td>DNA</td>
      <td>40</td>
      <td>47.5</td>
      <td>77.5</td>
      <td>2.59</td>
    </tr>
    <tr>
      <th>10</th>
      <td>N2-LF</td>
      <td>GGG GGC AAA TTG TGC AAT TTG</td>
      <td>CAA ATT GCA CAA TTT GCC CCC</td>
      <td>0.4</td>
      <td>50</td>
      <td>8</td>
      <td>1.4</td>
      <td>DNA</td>
      <td>21</td>
      <td>47.6</td>
      <td>66.1</td>
      <td>4.85</td>
    </tr>
    <tr>
      <th>11</th>
      <td>N2-LB</td>
      <td>CTT CGG GAA CGT GGT TGA CC</td>
      <td>GGT CAA CCA CGT TCC CGA AG</td>
      <td>0.4</td>
      <td>50</td>
      <td>8</td>
      <td>1.4</td>
      <td>DNA</td>
      <td>20</td>
      <td>60.0</td>
      <td>67.3</td>
      <td>5.37</td>
    </tr>
    <tr>
      <th>12</th>
      <td>ACTB-F3</td>
      <td>AGT ACC CCA TCG AGC ACG</td>
      <td>CGT GCT CGA TGG GGT ACT</td>
      <td>0.2</td>
      <td>50</td>
      <td>8</td>
      <td>1.4</td>
      <td>DNA</td>
      <td>18</td>
      <td>61.1</td>
      <td>64.8</td>
      <td>5.73</td>
    </tr>
    <tr>
      <th>13</th>
      <td>ACTB-B3</td>
      <td>AGC CTG GAT AGC AAC GTA CA</td>
      <td>TGT ACG TTG CTA TCC AGG CT</td>
      <td>0.2</td>
      <td>50</td>
      <td>8</td>
      <td>1.4</td>
      <td>DNA</td>
      <td>20</td>
      <td>50.0</td>
      <td>64.9</td>
      <td>4.93</td>
    </tr>
    <tr>
      <th>14</th>
      <td>ACTB-FIP</td>
      <td>GAG CCA CAC GCA GCT CAT TGT ATC ACC AAC TGG GA...</td>
      <td>TGT CGT CCC AGT TGG TGA TAC AAT GAG CTG CGT GT...</td>
      <td>1.6</td>
      <td>50</td>
      <td>8</td>
      <td>1.4</td>
      <td>DNA</td>
      <td>40</td>
      <td>55.0</td>
      <td>79.0</td>
      <td>2.59</td>
    </tr>
    <tr>
      <th>15</th>
      <td>ACTB-BIP</td>
      <td>CTG AAC CCC AAG GCC AAC CGG CTG GGG TGT TGA AG...</td>
      <td>GAC CTT CAA CAC CCC AGC CGG TTG GCC TTG GGG TT...</td>
      <td>1.6</td>
      <td>50</td>
      <td>8</td>
      <td>1.4</td>
      <td>DNA</td>
      <td>38</td>
      <td>63.2</td>
      <td>80.9</td>
      <td>2.79</td>
    </tr>
    <tr>
      <th>16</th>
      <td>ACTB-LF</td>
      <td>TGT GGT GCC AGA TTT TCT CCA</td>
      <td>TGG AGA AAA TCT GGC ACC ACA</td>
      <td>0.4</td>
      <td>50</td>
      <td>8</td>
      <td>1.4</td>
      <td>DNA</td>
      <td>21</td>
      <td>47.6</td>
      <td>66.9</td>
      <td>5.19</td>
    </tr>
    <tr>
      <th>17</th>
      <td>ACTB-LB</td>
      <td>CGA GAA GAT GAC CCA GAT CAT GT</td>
      <td>ACA TGA TCT GGG TCA TCT TCT CG</td>
      <td>0.4</td>
      <td>50</td>
      <td>8</td>
      <td>1.4</td>
      <td>DNA</td>
      <td>23</td>
      <td>47.8</td>
      <td>66.1</td>
      <td>4.28</td>
    </tr>
  </tbody>
</table>
</div>




```python
!mkdir "Output Data GENE Analysis"
```


```python
df1.to_csv('Output Data GENE Analysis/GENE_Analysis.csv', index=False)
```


```python
%%time

Primer_Name = []
Sequence = []
Thermo = []
DeltaS = []
DeltaG = []
DeltaH = []

count = 0

for i in tqdm(list(GENE_df.Sequence)):
    data = '{  "Sequence": "' + i + '",  "NaConc": 50,  "FoldingTemp": 37,\
    "MgConc": 8, "NucleotideType": "DNA" }'
    response = requests.post(
        'https://www.idtdna.com/Restapi/v1/OligoAnalyzer/Hairpin',
        headers=headers,
        data=data)
    json_data = json.loads(response.text)
    Primer_Name.append(Sequence_dict[json_data[0]['sequence']])
    Sequence.append(json_data[0]['sequence'])
    Thermo.append(json_data[0]['thermo'])
    DeltaS.append(json_data[0]['deltaS'])
    DeltaG.append(json_data[0]['deltaG'])
    DeltaH.append(json_data[0]['deltaH'])

df2 = pd.DataFrame({
    "Primer Name": Primer_Name,
    "Sequence": Sequence,
    "Thermo": Thermo,
    "DeltaG": DeltaG,
    "DeltaS": DeltaS,
    "DeltaH": DeltaH,
})


```


      0%|          | 0/18 [00:00<?, ?it/s]


    CPU times: user 412 ms, sys: 36.9 ms, total: 449 ms
    Wall time: 17.6 s
    

# **<span style="color:#4086f7; font-size:22px">If the highest hairpin Tm is at or above your annealing temperature, that hairpin is likely to impede hybridization</span>**


```python
df2
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Primer Name</th>
      <th>Sequence</th>
      <th>Thermo</th>
      <th>DeltaG</th>
      <th>DeltaS</th>
      <th>DeltaH</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Gene-E1-F3</td>
      <td>TGAGTACGAACTTATGTACTCAT</td>
      <td>45.9</td>
      <td>-3.59</td>
      <td>-172.08</td>
      <td>-54.9</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Gene-E1-B3</td>
      <td>TTCAGATTTTTAACACGAGAGT</td>
      <td>30.9</td>
      <td>-0.23</td>
      <td>-38.81</td>
      <td>-11.8</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Gene-E1-FIP</td>
      <td>ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG</td>
      <td>42.6</td>
      <td>-1.77</td>
      <td>-100.40</td>
      <td>-31.7</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Gene-E1-BIP</td>
      <td>TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT</td>
      <td>49.1</td>
      <td>-3.46</td>
      <td>-143.67</td>
      <td>-46.3</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Gene-E1-LF</td>
      <td>CGCTATTAACTATTAACG</td>
      <td>15.7</td>
      <td>0.67</td>
      <td>-72.02</td>
      <td>-20.8</td>
    </tr>
    <tr>
      <th>5</th>
      <td>Gene-E1-LB</td>
      <td>GCGCTTCGATTGTGTGCGT</td>
      <td>39.4</td>
      <td>-1.62</td>
      <td>-112.64</td>
      <td>-35.2</td>
    </tr>
    <tr>
      <th>6</th>
      <td>N2-F3</td>
      <td>ACCAGGAACTAATCAGACAAG</td>
      <td>40.8</td>
      <td>-0.39</td>
      <td>-24.52</td>
      <td>-7.7</td>
    </tr>
    <tr>
      <th>7</th>
      <td>N2-B3</td>
      <td>GACTTGATCTTTGAAATTTGGATCT</td>
      <td>29.8</td>
      <td>-0.63</td>
      <td>-130.72</td>
      <td>-39.6</td>
    </tr>
    <tr>
      <th>8</th>
      <td>N2-FIP</td>
      <td>TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC</td>
      <td>41.5</td>
      <td>-3.16</td>
      <td>-191.67</td>
      <td>-60.3</td>
    </tr>
    <tr>
      <th>9</th>
      <td>N2-BIP</td>
      <td>CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA</td>
      <td>34.1</td>
      <td>-1.08</td>
      <td>-118.45</td>
      <td>-36.4</td>
    </tr>
    <tr>
      <th>10</th>
      <td>N2-LF</td>
      <td>GGGGGCAAATTGTGCAATTTG</td>
      <td>43.1</td>
      <td>-2.88</td>
      <td>-159.70</td>
      <td>-50.5</td>
    </tr>
    <tr>
      <th>11</th>
      <td>N2-LB</td>
      <td>CTTCGGGAACGTGGTTGACC</td>
      <td>44.2</td>
      <td>-1.46</td>
      <td>-75.94</td>
      <td>-24.1</td>
    </tr>
    <tr>
      <th>12</th>
      <td>ACTB-F3</td>
      <td>AGTACCCCATCGAGCACG</td>
      <td>26.0</td>
      <td>-0.07</td>
      <td>-68.86</td>
      <td>-20.6</td>
    </tr>
    <tr>
      <th>13</th>
      <td>ACTB-B3</td>
      <td>AGCCTGGATAGCAACGTACA</td>
      <td>23.6</td>
      <td>0.10</td>
      <td>-73.46</td>
      <td>-21.8</td>
    </tr>
    <tr>
      <th>14</th>
      <td>ACTB-FIP</td>
      <td>GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA</td>
      <td>41.4</td>
      <td>-2.79</td>
      <td>-170.07</td>
      <td>-53.5</td>
    </tr>
    <tr>
      <th>15</th>
      <td>ACTB-BIP</td>
      <td>CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC</td>
      <td>63.8</td>
      <td>-6.65</td>
      <td>-171.22</td>
      <td>-57.7</td>
    </tr>
    <tr>
      <th>16</th>
      <td>ACTB-LF</td>
      <td>TGTGGTGCCAGATTTTCTCCA</td>
      <td>17.2</td>
      <td>0.67</td>
      <td>-85.77</td>
      <td>-24.9</td>
    </tr>
    <tr>
      <th>17</th>
      <td>ACTB-LB</td>
      <td>CGAGAAGATGACCCAGATCATGT</td>
      <td>34.4</td>
      <td>-0.82</td>
      <td>-86.48</td>
      <td>-26.6</td>
    </tr>
  </tbody>
</table>
</div>




```python
df2.to_csv('Output Data GENE Analysis/GENE_HairPins.csv', index=False)
```


```python
%%time

Primer_Name = []
Sequence_Bonds = []
Sequence_Sequences = []
Sequence_DeltaG = []
Sequence_BasePairs = []
Sequence_Dimer = []
Sequence_SequencePair = []
count = 1
for i in tqdm(list(GENE_df.Sequence)):
    #     print(str(i))
    params = {
        'primary': str(i),
    }
    temp = 1
    response = requests.post(
        'https://www.idtdna.com/Restapi/v1/OligoAnalyzer/SelfDimer',
        params=params,
        headers=headers)
    json_data_Sequence = json.loads(response.text)
    Primer_Name.append(Sequence_dict[i])
    Sequence_Sequences.append(i)
    Sequence_DeltaG.append(json_data_Sequence[0]['DeltaG'])
    Sequence_BasePairs.append(json_data_Sequence[0]['BasePairs'])
    Sequence_Dimer.append(json_data_Sequence[0]['Dimer'])
    Sequence_Bonds.append(json_data_Sequence[0]['Bonds'])
    Sequence_SequencePair.append(temp)
    temp = temp + 1
    Sequence_Sequences.append(i)
    Primer_Name.append(Sequence_dict[i])
    Sequence_DeltaG.append(json_data_Sequence[1]['DeltaG'])
    Sequence_BasePairs.append(json_data_Sequence[1]['BasePairs'])
    Sequence_Dimer.append(json_data_Sequence[1]['Dimer'])
    Sequence_Bonds.append(json_data_Sequence[1]['Bonds'])
    Sequence_SequencePair.append(temp)
    temp = 1



df3=pd.DataFrame({"Primer Name":Primer_Name,'Sequence Pair Number':Sequence_SequencePair,'Sequence':Sequence_Sequences,'DeltaG':Sequence_DeltaG,\
                  'BasePairs':Sequence_BasePairs,'Dimer':Sequence_Dimer,'Bonds':Sequence_Bonds,})
df3
```


      0%|          | 0/18 [00:00<?, ?it/s]


    CPU times: user 403 ms, sys: 28.4 ms, total: 431 ms
    Wall time: 19.1 s
    




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Primer Name</th>
      <th>Sequence Pair Number</th>
      <th>Sequence</th>
      <th>DeltaG</th>
      <th>BasePairs</th>
      <th>Dimer</th>
      <th>Bonds</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Gene-E1-F3</td>
      <td>1</td>
      <td>TGAGTACGAACTTATGTACTCAT</td>
      <td>-8.77</td>
      <td>7</td>
      <td>None</td>
      <td>[2, 2, 2, 2, 2, 2, 2, 0, 0, 1, 0, 0, 1, 0, 0, ...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Gene-E1-F3</td>
      <td>2</td>
      <td>TGAGTACGAACTTATGTACTCAT</td>
      <td>-3.65</td>
      <td>4</td>
      <td>None</td>
      <td>[1, 0, 0, 2, 2, 2, 2, 0, 0, 1, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Gene-E1-B3</td>
      <td>1</td>
      <td>TTCAGATTTTTAACACGAGAGT</td>
      <td>-4.85</td>
      <td>4</td>
      <td>None</td>
      <td>[0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 2, 2, 2, 0, 1, ...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Gene-E1-B3</td>
      <td>2</td>
      <td>TTCAGATTTTTAACACGAGAGT</td>
      <td>-3.61</td>
      <td>2</td>
      <td>None</td>
      <td>[0, 0, 0, 1, 0, 2, 2, 0, 1, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Gene-E1-FIP</td>
      <td>1</td>
      <td>ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG</td>
      <td>-9.08</td>
      <td>5</td>
      <td>None</td>
      <td>[0, 0, 0, 0, 2, 2, 2, 2, 2, 0, 0, 1, 1, 0, 0, ...</td>
    </tr>
    <tr>
      <th>5</th>
      <td>Gene-E1-FIP</td>
      <td>2</td>
      <td>ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG</td>
      <td>-8.47</td>
      <td>5</td>
      <td>None</td>
      <td>[0, 0, 0, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>6</th>
      <td>Gene-E1-BIP</td>
      <td>1</td>
      <td>TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT</td>
      <td>-8.64</td>
      <td>6</td>
      <td>None</td>
      <td>[1, 0, 2, 2, 2, 2, 2, 2, 0, 0, 0, 1, 1, 1, 1, ...</td>
    </tr>
    <tr>
      <th>7</th>
      <td>Gene-E1-BIP</td>
      <td>2</td>
      <td>TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT</td>
      <td>-6.30</td>
      <td>4</td>
      <td>None</td>
      <td>[2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>8</th>
      <td>Gene-E1-LF</td>
      <td>1</td>
      <td>CGCTATTAACTATTAACG</td>
      <td>-4.85</td>
      <td>4</td>
      <td>None</td>
      <td>[0, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>9</th>
      <td>Gene-E1-LF</td>
      <td>2</td>
      <td>CGCTATTAACTATTAACG</td>
      <td>-4.85</td>
      <td>4</td>
      <td>None</td>
      <td>[0, 0, 2, 2, 2, 2, 0, 0, 0, 1, 1, 1, 1, 0, 0, ...</td>
    </tr>
    <tr>
      <th>10</th>
      <td>Gene-E1-LB</td>
      <td>1</td>
      <td>GCGCTTCGATTGTGTGCGT</td>
      <td>-9.89</td>
      <td>4</td>
      <td>None</td>
      <td>[2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>11</th>
      <td>Gene-E1-LB</td>
      <td>2</td>
      <td>GCGCTTCGATTGTGTGCGT</td>
      <td>-6.76</td>
      <td>4</td>
      <td>None</td>
      <td>[0, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>12</th>
      <td>N2-F3</td>
      <td>1</td>
      <td>ACCAGGAACTAATCAGACAAG</td>
      <td>-3.07</td>
      <td>2</td>
      <td>None</td>
      <td>[0, 2, 2, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>13</th>
      <td>N2-F3</td>
      <td>2</td>
      <td>ACCAGGAACTAATCAGACAAG</td>
      <td>-1.60</td>
      <td>2</td>
      <td>None</td>
      <td>[1, 0, 0, 2, 2, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, ...</td>
    </tr>
    <tr>
      <th>14</th>
      <td>N2-B3</td>
      <td>1</td>
      <td>GACTTGATCTTTGAAATTTGGATCT</td>
      <td>-9.25</td>
      <td>6</td>
      <td>None</td>
      <td>[0, 0, 0, 1, 0, 0, 2, 2, 2, 2, 2, 2, 0, 0, 1, ...</td>
    </tr>
    <tr>
      <th>15</th>
      <td>N2-B3</td>
      <td>2</td>
      <td>GACTTGATCTTTGAAATTTGGATCT</td>
      <td>-4.62</td>
      <td>4</td>
      <td>None</td>
      <td>[0, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>16</th>
      <td>N2-FIP</td>
      <td>1</td>
      <td>TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC</td>
      <td>-10.20</td>
      <td>5</td>
      <td>None</td>
      <td>[2, 2, 2, 2, 2, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, ...</td>
    </tr>
    <tr>
      <th>17</th>
      <td>N2-FIP</td>
      <td>2</td>
      <td>TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC</td>
      <td>-9.28</td>
      <td>4</td>
      <td>None</td>
      <td>[2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>18</th>
      <td>N2-BIP</td>
      <td>1</td>
      <td>CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA</td>
      <td>-5.38</td>
      <td>4</td>
      <td>None</td>
      <td>[0, 0, 0, 1, 0, 1, 0, 0, 2, 2, 2, 2, 0, 0, 1, ...</td>
    </tr>
    <tr>
      <th>19</th>
      <td>N2-BIP</td>
      <td>2</td>
      <td>CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA</td>
      <td>-5.37</td>
      <td>4</td>
      <td>None</td>
      <td>[1, 0, 0, 2, 2, 2, 2, 0, 0, 1, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>20</th>
      <td>N2-LF</td>
      <td>1</td>
      <td>GGGGGCAAATTGTGCAATTTG</td>
      <td>-11.22</td>
      <td>7</td>
      <td>None</td>
      <td>[2, 2, 2, 2, 2, 2, 2, 0, 0, 1, 1, 1, 1, 1, 1, ...</td>
    </tr>
    <tr>
      <th>21</th>
      <td>N2-LF</td>
      <td>2</td>
      <td>GGGGGCAAATTGTGCAATTTG</td>
      <td>-7.05</td>
      <td>4</td>
      <td>None</td>
      <td>[0, 1, 0, 0, 0, 2, 2, 2, 2, 0, 0, 0, 1, 0, 0, ...</td>
    </tr>
    <tr>
      <th>22</th>
      <td>N2-LB</td>
      <td>1</td>
      <td>CTTCGGGAACGTGGTTGACC</td>
      <td>-6.30</td>
      <td>4</td>
      <td>None</td>
      <td>[0, 0, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>23</th>
      <td>N2-LB</td>
      <td>2</td>
      <td>CTTCGGGAACGTGGTTGACC</td>
      <td>-4.41</td>
      <td>3</td>
      <td>None</td>
      <td>[2, 2, 2, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>24</th>
      <td>ACTB-F3</td>
      <td>1</td>
      <td>AGTACCCCATCGAGCACG</td>
      <td>-6.76</td>
      <td>4</td>
      <td>None</td>
      <td>[1, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0, 0, 0, 1, 0, ...</td>
    </tr>
    <tr>
      <th>25</th>
      <td>ACTB-F3</td>
      <td>2</td>
      <td>AGTACCCCATCGAGCACG</td>
      <td>-3.65</td>
      <td>4</td>
      <td>None</td>
      <td>[0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>26</th>
      <td>ACTB-B3</td>
      <td>1</td>
      <td>AGCCTGGATAGCAACGTACA</td>
      <td>-6.30</td>
      <td>4</td>
      <td>None</td>
      <td>[0, 0, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>27</th>
      <td>ACTB-B3</td>
      <td>2</td>
      <td>AGCCTGGATAGCAACGTACA</td>
      <td>-3.65</td>
      <td>4</td>
      <td>None</td>
      <td>[0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>28</th>
      <td>ACTB-FIP</td>
      <td>1</td>
      <td>GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA</td>
      <td>-6.34</td>
      <td>4</td>
      <td>None</td>
      <td>[1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2, 2, 2, ...</td>
    </tr>
    <tr>
      <th>29</th>
      <td>ACTB-FIP</td>
      <td>2</td>
      <td>GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA</td>
      <td>-6.31</td>
      <td>4</td>
      <td>None</td>
      <td>[2, 2, 2, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, ...</td>
    </tr>
    <tr>
      <th>30</th>
      <td>ACTB-BIP</td>
      <td>1</td>
      <td>CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC</td>
      <td>-12.50</td>
      <td>6</td>
      <td>None</td>
      <td>[1, 0, 0, 0, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>31</th>
      <td>ACTB-BIP</td>
      <td>2</td>
      <td>CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC</td>
      <td>-9.75</td>
      <td>4</td>
      <td>None</td>
      <td>[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, ...</td>
    </tr>
    <tr>
      <th>32</th>
      <td>ACTB-LF</td>
      <td>1</td>
      <td>TGTGGTGCCAGATTTTCTCCA</td>
      <td>-5.02</td>
      <td>3</td>
      <td>None</td>
      <td>[1, 0, 2, 2, 2, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>33</th>
      <td>ACTB-LF</td>
      <td>2</td>
      <td>TGTGGTGCCAGATTTTCTCCA</td>
      <td>-5.02</td>
      <td>3</td>
      <td>None</td>
      <td>[2, 2, 2, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, ...</td>
    </tr>
    <tr>
      <th>34</th>
      <td>ACTB-LB</td>
      <td>1</td>
      <td>CGAGAAGATGACCCAGATCATGT</td>
      <td>-5.38</td>
      <td>4</td>
      <td>None</td>
      <td>[0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>35</th>
      <td>ACTB-LB</td>
      <td>2</td>
      <td>CGAGAAGATGACCCAGATCATGT</td>
      <td>-5.00</td>
      <td>4</td>
      <td>None</td>
      <td>[1, 0, 2, 2, 2, 2, 0, 1, 0, 0, 1, 0, 1, 1, 1, ...</td>
    </tr>
  </tbody>
</table>
</div>




```python
df3.to_csv('Output Data GENE Analysis/GENE_SelfDimers.csv', index=False)
```


```python
GENE_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Primer name</th>
      <th>Sequence</th>
      <th>Synthesis scale</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Gene-E1-F3</td>
      <td>TGAGTACGAACTTATGTACTCAT</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Gene-E1-B3</td>
      <td>TTCAGATTTTTAACACGAGAGT</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Gene-E1-FIP</td>
      <td>ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Gene-E1-BIP</td>
      <td>TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Gene-E1-LF</td>
      <td>CGCTATTAACTATTAACG</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>5</th>
      <td>Gene-E1-LB</td>
      <td>GCGCTTCGATTGTGTGCGT</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>6</th>
      <td>N2-F3</td>
      <td>ACCAGGAACTAATCAGACAAG</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>7</th>
      <td>N2-B3</td>
      <td>GACTTGATCTTTGAAATTTGGATCT</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>8</th>
      <td>N2-FIP</td>
      <td>TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>9</th>
      <td>N2-BIP</td>
      <td>CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>10</th>
      <td>N2-LF</td>
      <td>GGGGGCAAATTGTGCAATTTG</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>11</th>
      <td>N2-LB</td>
      <td>CTTCGGGAACGTGGTTGACC</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>12</th>
      <td>ACTB-F3</td>
      <td>AGTACCCCATCGAGCACG</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>13</th>
      <td>ACTB-B3</td>
      <td>AGCCTGGATAGCAACGTACA</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>14</th>
      <td>ACTB-FIP</td>
      <td>GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>15</th>
      <td>ACTB-BIP</td>
      <td>CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>16</th>
      <td>ACTB-LF</td>
      <td>TGTGGTGCCAGATTTTCTCCA</td>
      <td>10 nm</td>
    </tr>
    <tr>
      <th>17</th>
      <td>ACTB-LB</td>
      <td>CGAGAAGATGACCCAGATCATGT</td>
      <td>10 nm</td>
    </tr>
  </tbody>
</table>
</div>




```python
def SequencePairs(arr, n):
    a=[]
    for i in range(n):
        for j in range(n):
            a.append((arr[i],arr[j]))    
    return a
 
list_GENE_SEQUENCE=GENE_df.Sequence
n = len(list_GENE_SEQUENCE)
 
SequencePairs_list = SequencePairs(list_GENE_SEQUENCE, n)


for i in SequencePairs_list:
    if i[0]==i[1]:
        SequencePairs_list.remove(i)
        
        
for i in SequencePairs_list:
    if (i[0],i[1]) in SequencePairs_list and (i[1],i[0]) in SequencePairs_list:
        SequencePairs_list.remove((i[1],i[0]))
print('All the possible Primer Dimer possible Sets are \n')        
print(SequencePairs_list)
```

    All the possible Primer Dimer possible Sets are 
    
    [('TGAGTACGAACTTATGTACTCAT', 'TTCAGATTTTTAACACGAGAGT'), ('TGAGTACGAACTTATGTACTCAT', 'ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG'), ('TGAGTACGAACTTATGTACTCAT', 'TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT'), ('TGAGTACGAACTTATGTACTCAT', 'CGCTATTAACTATTAACG'), ('TGAGTACGAACTTATGTACTCAT', 'GCGCTTCGATTGTGTGCGT'), ('TGAGTACGAACTTATGTACTCAT', 'ACCAGGAACTAATCAGACAAG'), ('TGAGTACGAACTTATGTACTCAT', 'GACTTGATCTTTGAAATTTGGATCT'), ('TGAGTACGAACTTATGTACTCAT', 'TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC'), ('TGAGTACGAACTTATGTACTCAT', 'CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA'), ('TGAGTACGAACTTATGTACTCAT', 'GGGGGCAAATTGTGCAATTTG'), ('TGAGTACGAACTTATGTACTCAT', 'CTTCGGGAACGTGGTTGACC'), ('TGAGTACGAACTTATGTACTCAT', 'AGTACCCCATCGAGCACG'), ('TGAGTACGAACTTATGTACTCAT', 'AGCCTGGATAGCAACGTACA'), ('TGAGTACGAACTTATGTACTCAT', 'GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA'), ('TGAGTACGAACTTATGTACTCAT', 'CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC'), ('TGAGTACGAACTTATGTACTCAT', 'TGTGGTGCCAGATTTTCTCCA'), ('TGAGTACGAACTTATGTACTCAT', 'CGAGAAGATGACCCAGATCATGT'), ('TTCAGATTTTTAACACGAGAGT', 'ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG'), ('TTCAGATTTTTAACACGAGAGT', 'TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT'), ('TTCAGATTTTTAACACGAGAGT', 'CGCTATTAACTATTAACG'), ('TTCAGATTTTTAACACGAGAGT', 'GCGCTTCGATTGTGTGCGT'), ('TTCAGATTTTTAACACGAGAGT', 'ACCAGGAACTAATCAGACAAG'), ('TTCAGATTTTTAACACGAGAGT', 'GACTTGATCTTTGAAATTTGGATCT'), ('TTCAGATTTTTAACACGAGAGT', 'TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC'), ('TTCAGATTTTTAACACGAGAGT', 'CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA'), ('TTCAGATTTTTAACACGAGAGT', 'GGGGGCAAATTGTGCAATTTG'), ('TTCAGATTTTTAACACGAGAGT', 'CTTCGGGAACGTGGTTGACC'), ('TTCAGATTTTTAACACGAGAGT', 'AGTACCCCATCGAGCACG'), ('TTCAGATTTTTAACACGAGAGT', 'AGCCTGGATAGCAACGTACA'), ('TTCAGATTTTTAACACGAGAGT', 'GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA'), ('TTCAGATTTTTAACACGAGAGT', 'CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC'), ('TTCAGATTTTTAACACGAGAGT', 'TGTGGTGCCAGATTTTCTCCA'), ('TTCAGATTTTTAACACGAGAGT', 'CGAGAAGATGACCCAGATCATGT'), ('ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG', 'TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT'), ('ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG', 'CGCTATTAACTATTAACG'), ('ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG', 'GCGCTTCGATTGTGTGCGT'), ('ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG', 'ACCAGGAACTAATCAGACAAG'), ('ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG', 'GACTTGATCTTTGAAATTTGGATCT'), ('ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG', 'TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC'), ('ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG', 'CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA'), ('ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG', 'GGGGGCAAATTGTGCAATTTG'), ('ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG', 'CTTCGGGAACGTGGTTGACC'), ('ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG', 'AGTACCCCATCGAGCACG'), ('ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG', 'AGCCTGGATAGCAACGTACA'), ('ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG', 'GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA'), ('ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG', 'CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC'), ('ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG', 'TGTGGTGCCAGATTTTCTCCA'), ('ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG', 'CGAGAAGATGACCCAGATCATGT'), ('TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT', 'CGCTATTAACTATTAACG'), ('TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT', 'GCGCTTCGATTGTGTGCGT'), ('TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT', 'ACCAGGAACTAATCAGACAAG'), ('TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT', 'GACTTGATCTTTGAAATTTGGATCT'), ('TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT', 'TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC'), ('TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT', 'CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA'), ('TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT', 'GGGGGCAAATTGTGCAATTTG'), ('TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT', 'CTTCGGGAACGTGGTTGACC'), ('TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT', 'AGTACCCCATCGAGCACG'), ('TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT', 'AGCCTGGATAGCAACGTACA'), ('TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT', 'GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA'), ('TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT', 'CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC'), ('TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT', 'TGTGGTGCCAGATTTTCTCCA'), ('TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT', 'CGAGAAGATGACCCAGATCATGT'), ('CGCTATTAACTATTAACG', 'GCGCTTCGATTGTGTGCGT'), ('CGCTATTAACTATTAACG', 'ACCAGGAACTAATCAGACAAG'), ('CGCTATTAACTATTAACG', 'GACTTGATCTTTGAAATTTGGATCT'), ('CGCTATTAACTATTAACG', 'TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC'), ('CGCTATTAACTATTAACG', 'CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA'), ('CGCTATTAACTATTAACG', 'GGGGGCAAATTGTGCAATTTG'), ('CGCTATTAACTATTAACG', 'CTTCGGGAACGTGGTTGACC'), ('CGCTATTAACTATTAACG', 'AGTACCCCATCGAGCACG'), ('CGCTATTAACTATTAACG', 'AGCCTGGATAGCAACGTACA'), ('CGCTATTAACTATTAACG', 'GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA'), ('CGCTATTAACTATTAACG', 'CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC'), ('CGCTATTAACTATTAACG', 'TGTGGTGCCAGATTTTCTCCA'), ('CGCTATTAACTATTAACG', 'CGAGAAGATGACCCAGATCATGT'), ('GCGCTTCGATTGTGTGCGT', 'ACCAGGAACTAATCAGACAAG'), ('GCGCTTCGATTGTGTGCGT', 'GACTTGATCTTTGAAATTTGGATCT'), ('GCGCTTCGATTGTGTGCGT', 'TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC'), ('GCGCTTCGATTGTGTGCGT', 'CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA'), ('GCGCTTCGATTGTGTGCGT', 'GGGGGCAAATTGTGCAATTTG'), ('GCGCTTCGATTGTGTGCGT', 'CTTCGGGAACGTGGTTGACC'), ('GCGCTTCGATTGTGTGCGT', 'AGTACCCCATCGAGCACG'), ('GCGCTTCGATTGTGTGCGT', 'AGCCTGGATAGCAACGTACA'), ('GCGCTTCGATTGTGTGCGT', 'GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA'), ('GCGCTTCGATTGTGTGCGT', 'CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC'), ('GCGCTTCGATTGTGTGCGT', 'TGTGGTGCCAGATTTTCTCCA'), ('GCGCTTCGATTGTGTGCGT', 'CGAGAAGATGACCCAGATCATGT'), ('ACCAGGAACTAATCAGACAAG', 'GACTTGATCTTTGAAATTTGGATCT'), ('ACCAGGAACTAATCAGACAAG', 'TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC'), ('ACCAGGAACTAATCAGACAAG', 'CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA'), ('ACCAGGAACTAATCAGACAAG', 'GGGGGCAAATTGTGCAATTTG'), ('ACCAGGAACTAATCAGACAAG', 'CTTCGGGAACGTGGTTGACC'), ('ACCAGGAACTAATCAGACAAG', 'AGTACCCCATCGAGCACG'), ('ACCAGGAACTAATCAGACAAG', 'AGCCTGGATAGCAACGTACA'), ('ACCAGGAACTAATCAGACAAG', 'GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA'), ('ACCAGGAACTAATCAGACAAG', 'CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC'), ('ACCAGGAACTAATCAGACAAG', 'TGTGGTGCCAGATTTTCTCCA'), ('ACCAGGAACTAATCAGACAAG', 'CGAGAAGATGACCCAGATCATGT'), ('GACTTGATCTTTGAAATTTGGATCT', 'TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC'), ('GACTTGATCTTTGAAATTTGGATCT', 'CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA'), ('GACTTGATCTTTGAAATTTGGATCT', 'GGGGGCAAATTGTGCAATTTG'), ('GACTTGATCTTTGAAATTTGGATCT', 'CTTCGGGAACGTGGTTGACC'), ('GACTTGATCTTTGAAATTTGGATCT', 'AGTACCCCATCGAGCACG'), ('GACTTGATCTTTGAAATTTGGATCT', 'AGCCTGGATAGCAACGTACA'), ('GACTTGATCTTTGAAATTTGGATCT', 'GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA'), ('GACTTGATCTTTGAAATTTGGATCT', 'CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC'), ('GACTTGATCTTTGAAATTTGGATCT', 'TGTGGTGCCAGATTTTCTCCA'), ('GACTTGATCTTTGAAATTTGGATCT', 'CGAGAAGATGACCCAGATCATGT'), ('TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC', 'CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA'), ('TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC', 'GGGGGCAAATTGTGCAATTTG'), ('TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC', 'CTTCGGGAACGTGGTTGACC'), ('TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC', 'AGTACCCCATCGAGCACG'), ('TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC', 'AGCCTGGATAGCAACGTACA'), ('TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC', 'GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA'), ('TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC', 'CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC'), ('TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC', 'TGTGGTGCCAGATTTTCTCCA'), ('TTCCGAAGAACGCTGAAGCGGAACTGATTACAAACATTGGCC', 'CGAGAAGATGACCCAGATCATGT'), ('CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA', 'GGGGGCAAATTGTGCAATTTG'), ('CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA', 'CTTCGGGAACGTGGTTGACC'), ('CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA', 'AGTACCCCATCGAGCACG'), ('CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA', 'AGCCTGGATAGCAACGTACA'), ('CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA', 'GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA'), ('CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA', 'CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC'), ('CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA', 'TGTGGTGCCAGATTTTCTCCA'), ('CGCATTGGCATGGAAGTCACAATTTGATGGCACCTGTGTA', 'CGAGAAGATGACCCAGATCATGT'), ('GGGGGCAAATTGTGCAATTTG', 'CTTCGGGAACGTGGTTGACC'), ('GGGGGCAAATTGTGCAATTTG', 'AGTACCCCATCGAGCACG'), ('GGGGGCAAATTGTGCAATTTG', 'AGCCTGGATAGCAACGTACA'), ('GGGGGCAAATTGTGCAATTTG', 'GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA'), ('GGGGGCAAATTGTGCAATTTG', 'CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC'), ('GGGGGCAAATTGTGCAATTTG', 'TGTGGTGCCAGATTTTCTCCA'), ('GGGGGCAAATTGTGCAATTTG', 'CGAGAAGATGACCCAGATCATGT'), ('CTTCGGGAACGTGGTTGACC', 'AGTACCCCATCGAGCACG'), ('CTTCGGGAACGTGGTTGACC', 'AGCCTGGATAGCAACGTACA'), ('CTTCGGGAACGTGGTTGACC', 'GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA'), ('CTTCGGGAACGTGGTTGACC', 'CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC'), ('CTTCGGGAACGTGGTTGACC', 'TGTGGTGCCAGATTTTCTCCA'), ('CTTCGGGAACGTGGTTGACC', 'CGAGAAGATGACCCAGATCATGT'), ('AGTACCCCATCGAGCACG', 'AGCCTGGATAGCAACGTACA'), ('AGTACCCCATCGAGCACG', 'GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA'), ('AGTACCCCATCGAGCACG', 'CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC'), ('AGTACCCCATCGAGCACG', 'TGTGGTGCCAGATTTTCTCCA'), ('AGTACCCCATCGAGCACG', 'CGAGAAGATGACCCAGATCATGT'), ('AGCCTGGATAGCAACGTACA', 'GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA'), ('AGCCTGGATAGCAACGTACA', 'CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC'), ('AGCCTGGATAGCAACGTACA', 'TGTGGTGCCAGATTTTCTCCA'), ('AGCCTGGATAGCAACGTACA', 'CGAGAAGATGACCCAGATCATGT'), ('GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA', 'CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC'), ('GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA', 'TGTGGTGCCAGATTTTCTCCA'), ('GAGCCACACGCAGCTCATTGTATCACCAACTGGGACGACA', 'CGAGAAGATGACCCAGATCATGT'), ('CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC', 'TGTGGTGCCAGATTTTCTCCA'), ('CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC', 'CGAGAAGATGACCCAGATCATGT'), ('TGTGGTGCCAGATTTTCTCCA', 'CGAGAAGATGACCCAGATCATGT')]
    


```python
%%time

Sequence_Name1 = []
Sequence_Name2 = []
Sequence_Bonds = []
Sequence_Sequences1 = []
Sequence_Sequences2 = []
Sequence_DeltaG = []
Sequence_BasePairs = []
Sequence_Dimer = []
Sequence_SequencePair = []
count = 1
for i in tqdm(SequencePairs_list):
    temp=0
    params = {
        'primary': i[0],
        'secondary': i[1],}
    response = requests.post('https://www.idtdna.com/Restapi/v1/OligoAnalyzer/HeteroDimer', params=params, headers=headers)

    json_data_Sequence = json.loads(response.text)
    
    Sequence_Name1.append(Sequence_dict[i[0]])
    Sequence_Name2.append(Sequence_dict[i[1]])
    Sequence_Sequences1.append(i[0])
    Sequence_Sequences2.append(i[1])
    Sequence_DeltaG.append(json_data_Sequence[0]['DeltaG'])
    Sequence_BasePairs.append(json_data_Sequence[0]['BasePairs'])
    Sequence_Dimer.append(json_data_Sequence[0]['Dimer'])
    Sequence_Bonds.append(json_data_Sequence[0]['Bonds'])
    Sequence_SequencePair.append(temp)
    temp = temp + 1
    Sequence_Name1.append(Sequence_dict[i[0]])
    Sequence_Name2.append(Sequence_dict[i[1]])
    Sequence_Sequences1.append(i[0])
    Sequence_Sequences2.append(i[1])
    Sequence_DeltaG.append(json_data_Sequence[1]['DeltaG'])
    Sequence_BasePairs.append(json_data_Sequence[1]['BasePairs'])
    Sequence_Dimer.append(json_data_Sequence[1]['Dimer'])
    Sequence_Bonds.append(json_data_Sequence[1]['Bonds'])
    Sequence_SequencePair.append(temp)
    temp = 1
```


      0%|          | 0/153 [00:00<?, ?it/s]


    CPU times: user 3.21 s, sys: 225 ms, total: 3.44 s
    Wall time: 2min 42s
    

- **<span style="color:#4086f7; font-size:20px">The Delta G value of any heterodimers should be weaker (more positive) than -9.0 kcal/mole.</span>** 
- **<span style="color:#4086f7; font-size:20px">Positive numbers indicate that the actual secondary structure shown will not form at all.</span>**


```python
df4=pd.DataFrame({'Primary Sequence name':Sequence_Name1,'Secondary Sequence name':Sequence_Name2,\
                  'Sequence Pair Number':Sequence_SequencePair,'Primary Sequence':Sequence_Sequences1,\
                  'Secondary Sequence':Sequence_Sequences2,'DeltaG':Sequence_DeltaG,\
                  'BasePairs':Sequence_BasePairs,'Dimer':Sequence_Dimer,'Bonds':Sequence_Bonds,})
df4
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Primary Sequence name</th>
      <th>Secondary Sequence name</th>
      <th>Sequence Pair Number</th>
      <th>Primary Sequence</th>
      <th>Secondary Sequence</th>
      <th>DeltaG</th>
      <th>BasePairs</th>
      <th>Dimer</th>
      <th>Bonds</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Gene-E1-F3</td>
      <td>Gene-E1-B3</td>
      <td>0</td>
      <td>TGAGTACGAACTTATGTACTCAT</td>
      <td>TTCAGATTTTTAACACGAGAGT</td>
      <td>-4.52</td>
      <td>4</td>
      <td>None</td>
      <td>[2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Gene-E1-F3</td>
      <td>Gene-E1-B3</td>
      <td>1</td>
      <td>TGAGTACGAACTTATGTACTCAT</td>
      <td>TTCAGATTTTTAACACGAGAGT</td>
      <td>-3.61</td>
      <td>2</td>
      <td>None</td>
      <td>[0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 1, 0, 1, 0, 0, ...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Gene-E1-F3</td>
      <td>Gene-E1-FIP</td>
      <td>0</td>
      <td>TGAGTACGAACTTATGTACTCAT</td>
      <td>ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG</td>
      <td>-13.36</td>
      <td>8</td>
      <td>None</td>
      <td>[0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 0, 1, ...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Gene-E1-F3</td>
      <td>Gene-E1-FIP</td>
      <td>1</td>
      <td>TGAGTACGAACTTATGTACTCAT</td>
      <td>ACCACGAAAGCAAGAAAAAGAAGTTCGTTTCGGAAGAGACAG</td>
      <td>-7.13</td>
      <td>4</td>
      <td>None</td>
      <td>[0, 0, 0, 0, 1, 0, 2, 2, 2, 2, 0, 0, 0, 1, 0, ...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Gene-E1-F3</td>
      <td>Gene-E1-BIP</td>
      <td>0</td>
      <td>TGAGTACGAACTTATGTACTCAT</td>
      <td>TTGCTAGTTACACTAGCCATCCTTAGGTTTTACAAGACTCACGT</td>
      <td>-6.47</td>
      <td>5</td>
      <td>None</td>
      <td>[2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, ...</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>301</th>
      <td>ACTB-BIP</td>
      <td>ACTB-LF</td>
      <td>1</td>
      <td>CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC</td>
      <td>TGTGGTGCCAGATTTTCTCCA</td>
      <td>-6.21</td>
      <td>3</td>
      <td>None</td>
      <td>[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 2, 2, 0, ...</td>
    </tr>
    <tr>
      <th>302</th>
      <td>ACTB-BIP</td>
      <td>ACTB-LB</td>
      <td>0</td>
      <td>CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC</td>
      <td>CGAGAAGATGACCCAGATCATGT</td>
      <td>-9.69</td>
      <td>5</td>
      <td>None</td>
      <td>[0, 0, 1, 0, 0, 0, 0, 2, 2, 2, 2, 2, 0, 0, 0, ...</td>
    </tr>
    <tr>
      <th>303</th>
      <td>ACTB-BIP</td>
      <td>ACTB-LB</td>
      <td>1</td>
      <td>CTGAACCCCAAGGCCAACCGGCTGGGGTGTTGAAGGTC</td>
      <td>CGAGAAGATGACCCAGATCATGT</td>
      <td>-7.48</td>
      <td>4</td>
      <td>None</td>
      <td>[1, 0, 0, 0, 1, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0, ...</td>
    </tr>
    <tr>
      <th>304</th>
      <td>ACTB-LF</td>
      <td>ACTB-LB</td>
      <td>0</td>
      <td>TGTGGTGCCAGATTTTCTCCA</td>
      <td>CGAGAAGATGACCCAGATCATGT</td>
      <td>-6.69</td>
      <td>5</td>
      <td>None</td>
      <td>[1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 2, ...</td>
    </tr>
    <tr>
      <th>305</th>
      <td>ACTB-LF</td>
      <td>ACTB-LB</td>
      <td>1</td>
      <td>TGTGGTGCCAGATTTTCTCCA</td>
      <td>CGAGAAGATGACCCAGATCATGT</td>
      <td>-5.02</td>
      <td>3</td>
      <td>None</td>
      <td>[1, 0, 2, 2, 2, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, ...</td>
    </tr>
  </tbody>
</table>
<p>306 rows × 9 columns</p>
</div>




```python
df4.to_csv('Output Data GENE Analysis/GENE_Hetro_Dimers.csv', index=False)
```

# **<span style="color:#9f51e7;">References </span>**
- [How can I check my PCR primers using the OligoAnalyzer® program to ensure there are no significant primer design issues?](https://sg.idtdna.com/pages/support/faqs/how-can-i-check-my-pcr-primers-using-the-oligoanalyzer-program-to-ensure-there-are-no-significant-primer-design-issues-)

- [Mathematical model to reduce loop mediated isothermal amplification (LAMP) false-positive diagnosis](https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/full/10.1002/elps.201900167?casa_token=uh1wVkHdoH8AAAAA%3AObqNtNKeoTlabgQxXSncJFoB5qsBiIT-nOm1eErxLafmVXCjfxEduo_rhhHNAOk3PD8Vicl39lSj9A)

- [Impact  of  Primer-Dimers  and  Self-Amplifying  Hairpins  on  Reverse  Transcription  
Loop-Mediated Isothermal Amplification Detection of Viral RNA ](https://pubmed.ncbi.nlm.nih.gov/29620773/)

- [Thermodynamics and kinetics guided probe design for uniformly sensitive and specific DNA hybridization without optimization](https://www.nature.com/articles/s41467-019-12593-9)




 ### This Notebook is Created by [**Owais Ahmad**](https://www.linkedin.com/in/owaiskhan9654/) for SARS-CoV-2 Gene Sequence Analysis.
 ### This is for Medical Research Purpose only
 
 
 
- **Email owaiskhan9654@gmail.com**
- **Contact +919515884381**

<center><h1 style = "font-size:25px;font-family: Comic Sans MS"> Feel free to comment if you have any queries:)</h1></center>

