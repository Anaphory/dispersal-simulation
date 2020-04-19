def cultural_distance(c1, c2):
    return bin(c1 ^ c2).count("1")
def similar_culture(c1, c2):
    return bin(c1 ^ c2).count("1") < 6
