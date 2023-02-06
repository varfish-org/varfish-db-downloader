(NR == 1) {
    print $0;
    next;
}
(NR > 1) {
    print $0 | "LC_ALL=C sort -k1,1V -k2,2n -k3,3n"
}
