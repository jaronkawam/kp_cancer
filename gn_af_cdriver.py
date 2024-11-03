import time
import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By

# webdriver_path = '/usr/local/bin/'
driver = webdriver.Chrome()

# Open the website
varID = "4-83462758-T-C" #chrom, pos, ref, var
url = "https://gnomad.broadinstitute.org/variant/" + varID + "?dataset=gnomad_r4" # version 4.1, GRCh38
driver.get(url)

# Wait for the page to load
time.sleep(1) # You might need to adjust the time depending on your connection

# Find the table element (adjust the selector as necessary)
table = driver.find_element(By.XPATH, '//table[contains(@class, "Table__BaseTable-sc-7fgtt2-0 sc-cIShpX KcAIQ")]')

# Extract the table headers
headers = [header.text for header in table.find_elements(By.XPATH, './/thead/tr/th')]

# Extract the table rows
rows = []
for row in table.find_elements(By.XPATH, './/tbody/tr'):
    rows.append([row.text.split('\n')[0]] + [cell.text for cell in row.find_elements(By.XPATH, './/td')])

# Create a DataFrame from the extracted data
df = pd.DataFrame(rows, columns=headers)

print(df)

# Close the web driver
driver.quit()