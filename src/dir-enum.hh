#pragma once

enum class Direction : unsigned short { pos_x = 0,
                                        neg_x,
                                        pos_y,
                                        neg_y,
                                        pos_z,
                                        neg_z };


inline Direction to_direction(int i) {
  switch (i) {
    case 0:
      return Direction::pos_x;
    case 1:
      return Direction::neg_x;
    case 2:
      return Direction::pos_y;
    case 3:
      return Direction::neg_y;
    case 4:
      return Direction::pos_z;
    case 5:
      return Direction::neg_z;
  }
  assert(false);
}
