from reedsolo import RSCodec
from robust_soliton import PRNG


chunk_num = 1494
prng = PRNG(K=chunk_num, delta=0.05, c=0.1, np=False)
payload_2_chunks = []
chunks_2_payload = [None] * chunk_num
# Forms a bidirectional bipartite graph between chunks and payloads
file_chunks = [None] * chunk_num  # Index of chunk starts from 0
payloads = []


def dna_str_2_bytearray(dna_str: str) -> bytearray:
    # A=0, C=1, G=2, T=3
    dna_byte = bytearray()
    for i in range(0, len(dna_str), 4):
        val = 0
        for j in range(0, 4):
            val *= 4
            if dna_str[i + j] == 'A':
                val += 0
            elif dna_str[i + j] == 'C':
                val += 1
            elif dna_str[i + j] == 'G':
                val += 2
            else:
                val += 3
        dna_byte.append(val)
    return dna_byte


def bytes_xor(a: bytes, b: bytes) -> bytes:
    result = bytearray()
    for a_byte, b_byte in zip(a, b):
        result.append(a_byte ^ b_byte)
    return bytes(result)


def add_payload(chunks: set, payload: bytes):
    payload_idx = len(payloads)  # Index of payload starts from 0
    seen_chunks = set()
    for chunk in chunks:
        if file_chunks[chunk] is not None:
            file_chunk = file_chunks[chunk]
            payload = bytes_xor(payload, file_chunk)
            seen_chunks.add(chunk)
    for chunk in seen_chunks:
        chunks.remove(chunk)

    for chunk in chunks:
        if chunks_2_payload[chunk] is None:
            chunks_2_payload[chunk] = set()
        chunks_2_payload[chunk].add(payload_idx)
    payload_2_chunks.append(chunks)
    payloads.append(payload)

    propagate(payload_idx)


def propagate(payload_idx: int):
    if len(payload_2_chunks[payload_idx]) != 1:
        return
    
    corresponding_chunk_idx = payload_2_chunks[payload_idx].pop()
    file_chunks[corresponding_chunk_idx] = payloads[payload_idx]
    payloads[payload_idx] = bytes_xor(
        payloads[payload_idx],
        file_chunks[corresponding_chunk_idx]
    )  # This should make the payload zero
    
    associated_payloads = chunks_2_payload[corresponding_chunk_idx]
    for associated_payload in associated_payloads:
        # Skip the payload that is propagating
        if associated_payload == payload_idx:
            continue
        payload_2_chunks[associated_payload].remove(corresponding_chunk_idx)
        payloads[associated_payload] = bytes_xor(
            payloads[associated_payload],
            file_chunks[corresponding_chunk_idx]
        )

    for associated_payload in associated_payloads:
        propagate(associated_payload)


if __name__ == '__main__':
    with open('50-SF.txt', 'r') as f:
        lines = f.readlines()
        f.close()

    sequences = []
    for line in lines:
        line_strip = line.strip()
        if len(line_strip) == 100:
            if line_strip not in sequences:
                sequences.append(line_strip)
    encoded = [dna_str_2_bytearray(x) for x in sequences]

    rs = RSCodec(5)  # RS Code = 20 nt, 1 byte = 4 nt
    droplet_bytes = []
    for seq in encoded:
        try:
            corrected, _, __ = rs.decode(seq)
            corrected = bytes(corrected)
            if corrected not in droplet_bytes:
                droplet_bytes.append(corrected)
        except:
            pass

    required_droplets = 0
    for droplet_byte in droplet_bytes:
        required_droplets += 1
        seed_bytes = droplet_byte[:4]  # Just hardcode it
        seed = 0
        for seed_byte in seed_bytes:
            seed = (seed << 8) + seed_byte
        payload = droplet_byte[4:]
        _, chunks = prng.get_src_blocks_wrap(seed=seed)
        add_payload(set(chunks), payload)
        if None not in file_chunks:
            break

    with open('50-SF-decoded.jpg', 'wb') as f:
        for file_chunk in file_chunks:
            f.write(file_chunk)

    print(f'Full file is recovered after receiving {required_droplets} droplets.')
