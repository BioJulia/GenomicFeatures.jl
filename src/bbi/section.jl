# Data Section
# ============

# Each section is compressed separately.
struct Section
    chromid::UInt32
    chromstart::UInt32
    chromend::UInt32
    offset::UInt64
    size::UInt64
end
