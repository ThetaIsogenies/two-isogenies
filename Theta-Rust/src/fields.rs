#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

// ========================================================
// Definitions of base fields GF(p) = Z / pZ
// Constants defined are for the macro and generated from
// gen_fp.sage
// ========================================================
pub mod Fp254 {
    const N: usize = 4;
    const BITLEN: usize = 254;
    const MODULUS: [u64; N] = [
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xB52F88A2BBB638F2,
        0x300C882522D1C193,
    ];
    const HALF_MODULUS: [u64; N] = [
        0x0000000000000000,
        0x8000000000000000,
        0xDA97C4515DDB1C79,
        0x180644129168E0C9,
    ];
    const R_VAL: [u64; N] = [
        0x0000000000000005,
        0x0000000000000000,
        0x761254D25570E341,
        0x0FC1574651E7381D,
    ];
    const MINUS_R_VAL: [u64; N] = [
        0xFFFFFFFFFFFFFFFA,
        0xFFFFFFFFFFFFFFFF,
        0x3F1D33D0664555B1,
        0x204B30DED0EA8976,
    ];
    const DR_VAL: [u64; N] = [
        0x000000000000000A,
        0x0000000000000000,
        0xEC24A9A4AAE1C682,
        0x1F82AE8CA3CE703A,
    ];
    const TR_VAL: [u64; N] = [
        0x000000000000000F,
        0x0000000000000000,
        0x6236FE770052A9C3,
        0x2F4405D2F5B5A858,
    ];
    const QR_VAL: [u64; N] = [
        0x0000000000000015,
        0x0000000000000000,
        0x2319CAA69A0D5411,
        0x0EF8D4F424CB1EE2,
    ];
    const R2_VAL: [u64; N] = [
        0xE68D0176D1DE3F65,
        0xD5C29C23725CE9C4,
        0x333FEFEC3ADA8206,
        0x1CB17200F07FFDEA,
    ];
    const P0I: u64 = 1;
    const TFIXDIV_VAL: [u64; N] = [
        0x9DB4AE4A664B1867,
        0xF6DA2E7D8CD19F7B,
        0x5538A8E3E2FB8DD6,
        0x2B17BAA50A8EB9A9,
    ];
    const TDEC_VAL: [u64; N] = [
        0xD5C29C23725CE9C4,
        0xA60DA3D7E77CC6E5,
        0xBA5A0AE3EAA05CD2,
        0x2B45B97EB36F799D,
    ];
    const WIN_LEN: usize = 5;
    const SQRT_EH: [u8; 26] = [
        6, 15, 28, 24, 22, 27, 21, 8, 17, 24, 23, 20, 22, 19, 12, 16, 3, 13, 17, 20, 4, 8, 4, 3, 0,
        3,
    ];
    const SQRT_EL: usize = 25;
    const FOURTH_ROOT_EH: [u8; 26] = [
        19, 7, 14, 12, 27, 29, 10, 20, 8, 28, 11, 10, 27, 9, 6, 24, 17, 22, 8, 10, 2, 4, 18, 1, 16,
        1,
    ];
    const FOURTH_ROOT_EL: usize = 25;
    const P1: u64 = 3224510612;
    const P1DIV_M: u64 = 6123856568576173817;

    crate::finitefield::fp_gen::define_fp_core! {}

    #[cfg(test)]
    mod tests {
        crate::finitefield::fp_gen::define_fp_tests! {}
    }
}

pub mod Fp381 {
    const N: usize = 6;
    const BITLEN: usize = 381;
    const MODULUS: [u64; N] = [
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x3B1DC5C89B33FFFF,
        0x56456901A34100F2,
        0x14156D5C237E7FA2,
    ];
    const HALF_MODULUS: [u64; N] = [
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x1D8EE2E44D9A0000,
        0x2B22B480D1A08079,
        0x0A0AB6AE11BF3FD1,
    ];
    const R_VAL: [u64; N] = [
        0x000000000000000C,
        0x0000000000000000,
        0x0000000000000000,
        0x3A9ABA98B9900000,
        0xF4BF13EC58F3F4A5,
        0x0EFEDFAE56120463,
    ];
    const MINUS_R_VAL: [u64; N] = [
        0xFFFFFFFFFFFFFFF3,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x00830B2FE1A3FFFF,
        0x618655154A4D0C4D,
        0x05168DADCD6C7B3E,
    ];
    const DR_VAL: [u64; N] = [
        0x0000000000000019,
        0x0000000000000000,
        0x0000000000000000,
        0x3A17AF68D7EC0000,
        0x9338BED70EA6E858,
        0x09E8520088A58925,
    ];
    const TR_VAL: [u64; N] = [
        0x0000000000000026,
        0x0000000000000000,
        0x0000000000000000,
        0x3994A438F6480000,
        0x31B269C1C459DC0B,
        0x04D1C452BB390DE7,
    ];
    const QR_VAL: [u64; N] = [
        0x0000000000000032,
        0x0000000000000000,
        0x0000000000000000,
        0x742F5ED1AFD80000,
        0x26717DAE1D4DD0B0,
        0x13D0A401114B124B,
    ];
    const R2_VAL: [u64; N] = [
        0xBCEB199683E904B5,
        0xE30EDA89FF9E18B1,
        0x8855E34AEEB37EC9,
        0x79AE4E2B66A15CA7,
        0xE2D0EB9F2BFF4B66,
        0x0044A0A5FF0EFCBA,
    ];
    const P0I: u64 = 1;
    const TFIXDIV_VAL: [u64; N] = [
        0x3BDB13F19A3063FC,
        0xBAE0410E03C98849,
        0xC260B19F76D2E5EB,
        0x27B15CABF0D69F25,
        0x02795478145156CE,
        0x08C3645CD9D3A94E,
    ];
    const TDEC_VAL: [u64; N] = [
        0xE30EDA89FF9E18B1,
        0x8855E34AEEB37EC9,
        0x6C4B25C1F2655CA7,
        0x23846142DA93ED54,
        0xC01C3B2B3AB51708,
        0x0ED22DFCF43ACD23,
    ];
    const WIN_LEN: usize = 5;
    const SQRT_EH: [u8; 35] = [
        8, 19, 13, 2, 25, 5, 14, 7, 22, 3, 25, 3, 0, 1, 26, 8, 3, 16, 20, 21, 8, 22, 18, 8, 31, 7,
        31, 13, 4, 28, 10, 27, 10, 1, 10,
    ];
    const SQRT_EL: usize = 41;
    const FOURTH_ROOT_EH: [u8; 35] = [
        20, 25, 6, 17, 28, 2, 23, 3, 27, 17, 28, 1, 16, 0, 13, 20, 1, 8, 26, 10, 4, 11, 9, 20, 31,
        19, 31, 6, 2, 14, 21, 13, 21, 0, 5,
    ];
    const FOURTH_ROOT_EL: usize = 41;
    const P1: u64 = 2695588577;
    const P1DIV_M: u64 = 10945041894770733124;

    crate::finitefield::fp_gen::define_fp_core! {}

    #[cfg(test)]
    mod tests {
        crate::finitefield::fp_gen::define_fp_tests! {}
    }
}

pub mod FpFESTA {
    const N: usize = 21;
    const BITLEN: usize = 1293;
    const MODULUS: [u64; N] = [
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x43FFFFFFFFFFFFFF,
        0x9AC3245C4D7BE6B3,
        0x21D7DCD797059B7B,
        0x8F19A73E323F6569,
        0x841FED4773CFDB16,
        0x02979D50DD13D09A,
        0x01712922BAF59934,
        0xBD1C756E54F72C15,
        0xF6B3CF47C54370FE,
        0xCEC87BD4C1480F2B,
        0x11CF13E54B11406F,
        0x000000000000176C,
    ];
    const HALF_MODULUS: [u64; N] = [
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0xA200000000000000,
        0xCD61922E26BDF359,
        0x90EBEE6BCB82CDBD,
        0x478CD39F191FB2B4,
        0x420FF6A3B9E7ED8B,
        0x014BCEA86E89E84D,
        0x80B894915D7ACC9A,
        0x5E8E3AB72A7B960A,
        0xFB59E7A3E2A1B87F,
        0xE7643DEA60A40795,
        0x08E789F2A588A037,
        0x0000000000000BB6,
    ];
    const R_VAL: [u64; N] = [
        0x000AEE091BECF8C7,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x2400000000000000,
        0x78FB641CC3C350C6,
        0xA7023CB1847A0623,
        0x81CA1B3E72AA42FF,
        0x375C61B4D97A8C2B,
        0x36365530C09FB122,
        0xC85A4827D5B4CC8C,
        0x07F34CE9CA1014FC,
        0x5515312398D806A2,
        0x1051CF6B62E881A8,
        0xD49418553ECD1FE6,
        0x0000000000000A37,
    ];
    const MINUS_R_VAL: [u64; N] = [
        0xFFF511F6E4130738,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x1FFFFFFFFFFFFFFF,
        0x21C7C03F89B895ED,
        0x7AD5A026128B9558,
        0x0D4F8BFFBF952269,
        0x4CC38B929A554EEB,
        0xCC6148201C741F78,
        0x3916E0FAE540CCA7,
        0xB52928848AE71718,
        0xA19E9E242C6B6A5C,
        0xBE76AC695E5F8D83,
        0x3D3AFB900C442089,
        0x0000000000000D34,
    ];
    const DR_VAL: [u64; N] = [
        0x0015DC1237D9F18E,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x4800000000000000,
        0xF1F6C8398786A18C,
        0x4E04796308F40C46,
        0x0394367CE55485FF,
        0x6EB8C369B2F51857,
        0x6C6CAA61813F6244,
        0x90B4904FAB699918,
        0x0FE699D3942029F9,
        0xAA2A624731B00D44,
        0x20A39ED6C5D10350,
        0xA92830AA7D9A3FCC,
        0x000000000000146F,
    ];
    const TR_VAL: [u64; N] = [
        0x0020CA1B53C6EA56,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x2800000000000000,
        0xD02F07F9FDCE0B9F,
        0xD32ED93CF66876EE,
        0xF644AA7D25BF6395,
        0x21F537D7189FC96B,
        0xA00B624164CB42CC,
        0x579DAF54C628CC70,
        0x5ABD714F093912E1,
        0x088BC4230544A2E7,
        0x622CF26D677175CD,
        0x6BED351A71561F42,
        0x000000000000073B,
    ];
    const QR_VAL: [u64; N] = [
        0x002BB8246FB3E31D,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x4C00000000000000,
        0x492A6C16C1915C65,
        0x7A3115EE7AE27D12,
        0x780EC5BB9869A695,
        0x5951998BF21A5597,
        0xD641B772256AF3EE,
        0x1FF7F77C9BDD98FC,
        0x62B0BE38D34927DE,
        0x5DA0F5469E1CA989,
        0x727EC1D8CA59F775,
        0x40814D6FB0233F28,
        0x0000000000001173,
    ];
    const R2_VAL: [u64; N] = [
        0x9F42DE4F788045C2,
        0x1C2EB545BBE1E6CA,
        0x2C7698FD862CBB36,
        0x00FDEE42F690BB02,
        0x37B26D78444A90CA,
        0x95622499BE0F470F,
        0xC62D0EAF8B7DC9CC,
        0xF242919A892AA2DF,
        0x601AD6354D69F3CB,
        0x074F84AEB6FA3F2B,
        0xF14BEC8FE0ABEB2D,
        0xD8339C56ACAD0170,
        0xF06FF606EE809357,
        0x88D41B5189622D47,
        0xB6C213020828FB73,
        0x9F8125EC514BBAAC,
        0xAC872B8F8CFA8AD7,
        0x711147A40E5E32CA,
        0xBCF167C3A3CBB9FE,
        0x76F7DA793EAA61EF,
        0x0000000000000D90,
    ];
    const P0I: u64 = 1;
    const TFIXDIV_VAL: [u64; N] = [
        0x57CFDEF044149D6F,
        0xCE41D70241B51263,
        0xFB982E57F54C7E28,
        0xA4679409D1E483C7,
        0x011DC2D5A0BAADE8,
        0x8AAA62A138C96A66,
        0xD222A8115270BA21,
        0x06C5848819FD1F77,
        0x9BB2C052D2153FD6,
        0x066E9B66267C36CA,
        0x3D5097316A500C8B,
        0xC61ADE5C138AE9CA,
        0x3C7F105387EA0D08,
        0x9179D77EB3F9D1DA,
        0x94B6937BF7B9E519,
        0xDF141B79C2C7EFA0,
        0xE30A188AF2ECCFE1,
        0xF7E4F400401E6695,
        0x7C3B040888C99764,
        0xDCF405E68822CC43,
        0x00000000000005D2,
    ];
    const TDEC_VAL: [u64; N] = [
        0x1C2EB545BBE1E6CA,
        0x2C7698FD862CBB36,
        0x00FDEE42F690BB02,
        0x37B26D78444A90CA,
        0x95622499BE0F470F,
        0xC62D0EAF8B7DC9CC,
        0xF242919A892AA2DF,
        0x601AD6354D69F3CB,
        0x8F4F84AEB6FA3F2B,
        0x71325BC65B41105A,
        0x79727C59FFDA3258,
        0xFF3C908B8ED06518,
        0xE1D183D9987FBC43,
        0x43072C6DA68FE17C,
        0x8A0380573C4B2FAD,
        0x108DB7C73303ED54,
        0xD97B27C377C7B8EA,
        0x792E18D4CD60A8A6,
        0x628D06B8251703F4,
        0x3D43319001A8F373,
        0x0000000000000E92,
    ];
    const WIN_LEN: usize = 5;
    const SQRT_EH: [u8; 133] = [
        4, 26, 12, 13, 30, 29, 21, 9, 28, 2, 9, 6, 12, 13, 14, 15, 27, 12, 1, 14, 25, 11, 19, 27,
        23, 14, 8, 18, 22, 18, 29, 7, 18, 17, 15, 14, 26, 12, 28, 17, 22, 24, 22, 31, 28, 25, 29,
        8, 13, 31, 7, 8, 8, 13, 2, 26, 19, 8, 23, 1, 21, 14, 30, 18, 2, 0, 13, 18, 25, 26, 11, 23,
        2, 9, 10, 2, 23, 0, 20, 2, 12, 25, 29, 9, 5, 23, 21, 14, 28, 8, 15, 29, 15, 24, 13, 8, 5,
        30, 17, 30, 28, 25, 26, 30, 11, 25, 3, 16, 20, 0, 19, 26, 27, 3, 18, 29, 28, 23, 1, 8, 17,
        24, 18, 10, 30, 9, 28, 25, 17, 0, 27, 14, 1,
    ];
    const SQRT_EL: usize = 126;
    const FOURTH_ROOT_EH: [u8; 132] = [
        2, 13, 22, 6, 31, 30, 26, 4, 14, 17, 4, 3, 22, 6, 23, 23, 13, 22, 0, 23, 28, 21, 25, 29,
        11, 7, 4, 9, 11, 25, 30, 3, 25, 24, 7, 7, 13, 6, 30, 8, 11, 12, 27, 15, 30, 28, 14, 20, 22,
        31, 3, 4, 20, 6, 1, 29, 9, 20, 27, 16, 10, 7, 15, 9, 1, 16, 6, 25, 12, 29, 21, 11, 17, 4,
        5, 17, 11, 0, 10, 1, 22, 28, 30, 20, 18, 27, 10, 7, 14, 20, 23, 30, 7, 28, 6, 20, 2, 31, 8,
        15, 30, 12, 13, 31, 21, 28, 1, 8, 10, 16, 9, 29, 29, 1, 25, 14, 30, 27, 0, 20, 8, 12, 9, 5,
        31, 4, 30, 28, 8, 16, 13, 23,
    ];
    const FOURTH_ROOT_EL: usize = 126;
    const P1: u64 = 3143667320;
    const P1DIV_M: u64 = 6755719943463976019;

    crate::finitefield::fp_gen::define_fp_core! {}

    #[cfg(test)]
    mod tests {
        crate::finitefield::fp_gen::define_fp_tests! {}
    }
}

// ========================================================
// Definitions of extension fields above the base fields
// GF(p^2) with modulus x^2 + 1 = 0 (using p = 3 mod 4)
// ========================================================

pub mod Fp254Ext {
    use super::Fp254::Fp;
    const NQR_RE: Fp = Fp::new([
        0x19067C2EE042CAE9,
        0x7FAA1A2A81781083,
        0xF0EE8CA2EBCDDCE4,
        0x2561C37B60330EF1,
    ]);

    crate::finitefield::fp2_gen::define_fp2_core! {}
    #[cfg(test)]
    mod tests {
        crate::finitefield::fp2_gen::define_fp2_tests! {}
    }
}

pub mod Fp381Ext {
    use super::Fp381::Fp;
    const NQR_RE: Fp = Fp::new([
        0xA5A3D7A6B1EF921B,
        0xBE06313F3FB66DF8,
        0xCAE1514B24E62622,
        0xA85848716103C84B,
        0xFBC52D7FC56C0B0D,
        0x0703DA7C0679837C,
    ]);

    crate::finitefield::fp2_gen::define_fp2_core! {}
    #[cfg(test)]
    mod tests {
        crate::finitefield::fp2_gen::define_fp2_tests! {}
    }
}

pub mod FpFESTAExt {
    use super::FpFESTA::Fp;
    const NQR_RE: Fp = Fp::new([
        0x145E4F94B005484D,
        0x437E234B4E5CF5D6,
        0x5D29F56B6150B469,
        0xD1DCEFD038BFD33D,
        0xC3C772943D654936,
        0xCEAC95C36C49B427,
        0xE282B9FB7D9CBE13,
        0x0677E9713A992BF4,
        0xB7ECEF2567197DCC,
        0x8E9451C2845DFB09,
        0xA1229BB35705E1B9,
        0x38BC621169F0A3C7,
        0x215B414B1AB53288,
        0x1C3384A030CC972B,
        0xD7E9BB1787954CE7,
        0xCD570AD9A9444C6F,
        0x3EC28136E459A065,
        0x96D89B4337340AC7,
        0xA92B8E9DEBBBEB2B,
        0xA4D7F84010F736DF,
        0x0000000000001047,
    ]);

    crate::finitefield::fp2_gen::define_fp2_core! {}
    #[cfg(test)]
    mod tests {
        // crate::finitefield::fp2_gen::define_fp2_tests! {}
    }
}
