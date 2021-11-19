#pragma once

#include "Image/ImageAlgorithms.h"

#include <iostream>
#include <memory>
#include <vector>

class ImgSlots {
public:
    std::shared_ptr<baseImage> get(size_t slot, bool remove = false);
    void set(size_t slot, std::shared_ptr<baseImage> img);

private:
    std::vector<std::shared_ptr<baseImage>> m_slots;
};

std::ostream& operator<<(std::ostream& stream, const baseImage& rhs);
