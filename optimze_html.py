from bs4 import BeautifulSoup
'''Removes the first two <script> tags and their contents as well as the <html> 
  <head> and <body> tags in order to optimize the files into the CFI webpage'''


proteins = ['Spike','Envelope','Membrane','Nucleocapsid']

for protein in proteins:
  for i in range(1,38):
    filename = '08_table_{}_{}.html'.format(protein,i)
    html = open('output/epitsurve/' + filename)
    soup = BeautifulSoup(html,'html.parser')
    soup.select('script')[0].extract()
    soup.select('script')[0].extract()
    div = soup.select('div')[0]
    with open('output/optimized/' + filename, mode="w",  encoding="utf8") as code:
      code.write(str(div.prettify()))
