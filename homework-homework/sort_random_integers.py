import random

numbers = [random.randint(1, 750) for _ in range(20)]

print("random numbers:")
print(numbers)

for i in range(len(numbers)):
    for j in range(len(numbers) - 1 - i):
        if numbers[j] > numbers[j + 1]:
            numbers[j], numbers[j + 1] = numbers[j + 1], numbers[j]

print("\nSorted numbers:")
print(numbers)