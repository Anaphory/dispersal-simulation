pub fn smaller_of_two<T: PartialOrd>(a: T, b: T) -> T {
    if a < b {
        a
    } else {
        b
    }
}

pub fn greater_of_two<T: PartialOrd>(a: T, b: T) -> T {
    if b < a {
        a
    } else {
        b
    }
}

