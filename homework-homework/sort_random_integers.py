import random
numbers = [random.randint(1,750) for _ in range(20)]
print("Original list:")
print(numbers)
for i in range(len(numbers)):
    switched = False
    for j in range(len(numbers) - 1 - i):
        if numbers[j] > numbers[j + 1]:
            temp = numbers[j]
            numbers[j] = numbers[j + 1]
            numbers[j + 1] = temp

            switched = True

    if not switched:
        break
print("Sorted list:")
print(numbers)