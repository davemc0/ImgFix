#include "Helpers.h"

std::shared_ptr<baseImage> ImgSlots::get(size_t slot)
{
    std::cerr << "Pointing curImage at same image as slot " << slot << ".\n";
    if (m_slots.size() <= slot || slot < 0) throw DMcError("Slot out of range: " + std::to_string(m_slots.size()));
    if (!m_slots[slot]) throw DMcError("Slot empty");

    std::shared_ptr<baseImage> tmp = m_slots[slot];
    m_slots[slot] = nullptr;

    return tmp;
}

void ImgSlots::set(size_t slot, std::shared_ptr<baseImage> img)
{
    std::cerr << "Storing curImage in slot " << slot << ".\n";
    if (slot < 0) throw DMcError("Slot out of range: " + std::to_string(m_slots.size()));
    if (m_slots.size() <= slot) m_slots.resize(slot + 1);

    m_slots[slot] = img;
}

std::ostream& operator<<(std::ostream& stream, const baseImage& rhs)
{
    stream << " size=" << rhs.w_virtual() << "x" << rhs.h_virtual() << (rhs.is_integer_virtual() ? " INT" : " FLOAT") << " chan=" << rhs.chan_virtual()
           << " bytes_per_element=" << rhs.size_element_virtual();
    return stream;
}
