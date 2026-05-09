def read(file1):
    matrix = []
    with open(file1) as file: 
        lines = []
        for line in file:
            line = line.strip()
            if line != "":
                lines.append(line)
    n = int(lines[0])
    for i in range(1, len(lines)):
        matrix.append(list(map(float, lines[i].split())))
    return n, matrix

def average_distance(cluster1, cluster2, matrix):
    total = 0 
    count = 0 
    for i in cluster1:
        for j in cluster2:
            total += matrix[i - 1][j - 1]
            count += 1
    return total / count
        
def Hierarchal_clustering(n, matrix):
    clusters = []
    for i in range(1, n +1):
        clusters.append([i])
    finals = []
    while len(clusters) > 1:
        best_point = 0 
        best_start = 1 
        best_distance = average_distance(clusters[0], clusters[1], matrix)
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                distance = average_distance(clusters[i], clusters[j], matrix)
                if distance < best_distance:
                    best_distance = distance
                    best_point = i 
                    best_start = j 
        new_cluster = clusters[best_point] + clusters[best_start]
        finals.append(new_cluster)
        if best_point < best_start:
            clusters.pop(best_start)
            clusters.pop(best_point)
        else: 
            clusters.pop(best_point)
            clusters.pop(best_start)
        clusters.append(new_cluster)
    return finals
if __name__ == "__main__":
    n, matrix = read("rosalind_ba8e.txt")
    finals = Hierarchal_clustering(n, matrix)
    for cluster in finals:
        values = []
        for i in cluster:
            values.append(str(i))
        print(" ".join(values))
