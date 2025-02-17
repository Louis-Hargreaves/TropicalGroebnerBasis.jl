#A function that takes a string and returns the number of words in the string.

function count_words(str::String)
    return length(split(str))
end