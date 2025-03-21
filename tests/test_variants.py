import pytest
from whatshap.variants import AlignedRead, Read
from hypothesis import given, strategies as st


def create_aligned_read(name, ref_start, ref_end):
    return AlignedRead(
        read=Read(name, 60, 0, 0, ref_start, "", -1, -1),
        is_supplementary=False,
        is_reverse=False,
        reference_start=ref_start,
        reference_end=ref_end,
    )


segment_strategy = st.tuples(
    st.integers(min_value=0, max_value=2**15), st.integers(min_value=0, max_value=2**15)
).filter(lambda x: x[0] <= x[1])


@given(
    segment=segment_strategy,
)
def test_identity_distance(segment):
    start, end = segment
    read = create_aligned_read("read", start, end)
    assert read.distance(read) == 0


@given(
    segment_a=segment_strategy,
    segment_b=segment_strategy,
)
def test_pairwise_distances(segment_a, segment_b):
    start_a, end_a = segment_a
    start_b, end_b = segment_b
    read_a = create_aligned_read("read_a", start_a, end_a)
    read_b = create_aligned_read("read_b", start_b, end_b)
    assert read_a.distance(read_b) >= 0
    assert read_b.distance(read_a) == read_a.distance(read_b)


@pytest.fixture
def reads():
    read0 = create_aligned_read("read0", 100, 200)
    read1 = create_aligned_read("read1", 150, 250)
    read2 = create_aligned_read("read2", 300, 400)
    read3 = create_aligned_read("read3", 200, 250)
    read4 = create_aligned_read("read4", 110, 120)
    return read0, read1, read2, read3, read4


@pytest.mark.parametrize(
    "index_a, index_b, expected_distance",
    [
        (0, 1, 0),
        (0, 2, 100),
        (0, 3, 0),
        (0, 4, 0),
        (1, 2, 50),
        (1, 3, 0),
        (1, 4, 30),
        (2, 3, 50),
        (2, 4, 180),
        (3, 4, 80),
    ],
)
def test_distance(reads, index_a, index_b, expected_distance):
    assert reads[index_a].distance(reads[index_b]) == expected_distance
