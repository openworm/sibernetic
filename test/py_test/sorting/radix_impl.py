from typing import List


def sort(l: List[int]) -> List[int]:
    k = max(l)
    num_bits = len(bin(k)[2:])
    zero_bags = []
    one_bags = []
    result = l
    for i in range(num_bits):
        for r in result:
            if not(r & (1 << i)):
                zero_bags.append(r)
            else:
                one_bags.append(r)
        result = zero_bags
        result += one_bags
        zero_bags = []
        one_bags = []
    return result


if __name__ == '__main__':
    l = [2, 3, 4, 5, 2, 3, 10, 1, 0, 22]
    print(sort(l))
