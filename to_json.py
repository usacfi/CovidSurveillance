from bs4 import BeautifulSoup
import os

proteins = ['Spike','Envelope','Membrane','Nucleocapsid']

for protein in proteins:
  for i in range(1,38):
    # mutation profile
    filename = '07_mutation_profile_{}_{}.html'.format(protein,i)
    html = open('output/optimized/' + filename)
    soup = BeautifulSoup(html,'html.parser')
    output = soup.get_text()
    loc = output.splitlines()[10].find('"title": {"text": "Sars-Cov-2')
    data = output.splitlines()[9][:-1].lstrip()
    layout = output.splitlines()[10][loc:]
    config = output.splitlines()[11].lstrip()
    
    path = 'output/epitsurve/data/{}/{}/'.format(protein,i).lower()
    if not(os.path.exists('output/epitsurve/data/{}/'.format(protein))):
      os.mkdir('output/epitsurve/data/{}/'.format(protein).lower())
    if not(os.path.exists(path)):  
      os.mkdir(path)
    open(path + 'graph-data.json','w').write(str(data))
    open(path + 'graph-layout.json','w').write(str('{' + '\n"layout":\n' + '{' + layout + '\n"config":\n' + config + '}'))
       
    # table
    filename = '08_table_{}_{}.html'.format(protein,i)
    html = open('output/optimized/' + filename)
    soup = BeautifulSoup(html,'html.parser')
    output = soup.get_text()
    loc = output.splitlines()[10].find('"title": {"text": "Summary Table')
    data = output.splitlines()[9].lstrip()
    layout = output.splitlines()[10][loc:]
    config = output.splitlines()[11].lstrip()
    
    open(path + 'table.json','w').write(str('{"data":\n' + data + '\n"layout":\n' + '{' + layout + '\n"config":\n' + config + '}'))
    
    
