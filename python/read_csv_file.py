import csv

csv_file = "python/dragen_runtime_with_output_size.csv"

columns = []
rows = []

# reading csv file
with open(csv_file, 'r') as csv_file:
    csv_reader = csv.reader(csv_file)

    columns = next(csv_reader)

    for row in csv_reader:
        rows.append(row)

    print("Total number of rows: %d" % (csv_reader.line_num))

print('Column names are: ' + ', '.join(column for column in columns))

print('\nFirst 10 rows are:\n')
for row in rows[:10]:
    for col in row:
        print("%10s" % col, end=" "),
    print('\n')
    