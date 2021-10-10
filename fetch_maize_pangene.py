from bs4 import BeautifulSoup
import urllib.request
import urllib.parse
import sys

gene = sys.argv[1]
out = sys.argv[2]

def get_html(url):
	send_headers = {
		'User-Agent':'Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/47.0.2526.80 Safari/537.36',
		'Accept':'*/*',
		'Connection':'keep-alive',

	}

	req = urllib.request.Request(url,headers = send_headers)
	response = urllib.request.urlopen(req)
	html = response.read().decode('utf-8')

	return html

def analyse(html):
	soup = BeautifulSoup(html,'lxml')
	"""
	for i in soup.find_all('a'):
		print(i.string)
	for i in soup.find_all('b'):
		print(i.string)
	"""
	text = soup.get_text().split("\n")
	#print (text)

	jud = False
	content = []
	for i in text:
		i = "".join(i.split())
		if i == "":
			pass
		else:
			#pass
			#if "\\" in i:
			if jud == True:
				content.append("".join(i.split()))
			if "".join(i.split()) == "computedby":
				jud = True


	return content

html = get_html("https://ajax1.maizegdb.org/record_data/gene_pangenome_data.php?id=" + gene + "&type=related_genemodels_pangenome")
#print (html)
content = analyse(html)
step = 4
b = [content[i:i+step] for i in range(0, len(content), step)]
fout = open(out, 'w')
for i in b:
	fout.write("\t".join(i) + "\n")

