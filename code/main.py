__author__ = 'yun'

import csv

from exprsearch import ExprSearch


def main():
    with open("seeds_dataset.txt") as csv_file:
        data_reader = csv.reader(csv_file, delimiter='\t')
        dataset = []
        for row in data_reader:
            if "?" not in row:
                dataset.append(row)
        search = ExprSearch(dataset[0], 4)
        best = search.search(dataset[1:])
        best.plot_data()


if __name__ == "__main__":
    main()
