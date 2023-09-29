def analyze_scRNA(data):
    # Your code here
     count = data.count("A")
     return count

def load_data():
    # Your code here
    data = "ATCG" * 5
    return data

def main():
    # Your code here
    data = load_data()
    result = analyze_scRNA(data)
    print(result)

if __name__ == "__main__":
    main()
