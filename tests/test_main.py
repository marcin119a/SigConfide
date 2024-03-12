from sigconfide.modelselection.analyzer import cosmic_fit

if __name__ == "__main__":
    cosmic_fit('data/reduced_data.dat', '.', mutation_count=1000)