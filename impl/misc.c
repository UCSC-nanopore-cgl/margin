/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "margin.h"

char *getTimeDescriptorFromSeconds(int64_t seconds) {
    int64_t minutes = (int64_t) (seconds / 60);
    int64_t hours = (int64_t) (minutes / 60);
    char *timeDescriptor;

    if (hours > 0) {
        timeDescriptor = stString_print("%"PRId64"h %"PRId64"m", hours,
                minutes - (hours * 60));
    } else if (minutes > 0) {
        timeDescriptor = stString_print("%"PRId64"m %"PRId64"s", minutes,
                seconds - (minutes * 60));
    } else {
        timeDescriptor = stString_print("%"PRId64"s", seconds);
    }
    return timeDescriptor;
}