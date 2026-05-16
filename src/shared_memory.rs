#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(unused_imports)]

use crate::*;

#[doc = "Original class `SharedMemory` at STAR/source/SharedMemory.h:85."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SharedMemory {
    pub key: i32,
    pub counter_key: i32,
    pub unload_last: bool,
    pub shm_id: i32,
    pub shared_counter_id: i32,
    pub counter_mem_attached: bool,
    pub mapped: bool,
    pub length: usize,
    pub is_allocator: bool,
    pub needs_allocation: bool,
    pub exception: SharedMemoryException,
    pub shared_objects_use_count_value: i32,
    pub clean_count: u32,
}

#[doc = "Original class `SharedMemoryException` at STAR/source/SharedMemory.h:33."]
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SharedMemoryException {
    pub error_code: i32,
    pub error_detail: i32,
}

#[doc = "Original `SharedMemory::SharedMemory` at STAR/source/SharedMemory.cpp:22. Args: key: key_t, unloadLast: bool"]
pub fn sharedmemory_l22_sharedmemory_sharedmemory(
    key: i32,
    unload_last: bool,
) -> crate::shared_memory::SharedMemory {
    let mut shared_memory = crate::shared_memory::SharedMemory {
        key,
        counter_key: key + 1,
        unload_last,
        shm_id: -1,
        shared_counter_id: -1,
        counter_mem_attached: false,
        mapped: false,
        length: 0,
        is_allocator: false,
        needs_allocation: true,
        exception: crate::shared_memory::SharedMemoryException::default(),
        shared_objects_use_count_value: 1,
        clean_count: 0,
    };
    let _ = sharedmemory_l244_sharedmemory_ensurecounter(&mut shared_memory);
    let _ = sharedmemory_l134_sharedmemory_openifexists(&mut shared_memory);
    shared_memory
}

#[doc = "Original `SharedMemory::~SharedMemory` at STAR/source/SharedMemory.cpp:37. Args: "]
pub fn sharedmemory_l37_sharedmemory_sharedmemory(
    shared_memory: &mut crate::shared_memory::SharedMemory,
) -> Result<String, crate::shared_memory::SharedMemoryException> {
    let in_use = sharedmemory_l277_sharedmemory_sharedobjectsusecount(shared_memory)? - 1;
    if let Err(exc) = sharedmemory_l191_sharedmemory_close(shared_memory) {
        let _ = sharedmemory_l237_sharedmemory_clean(shared_memory);
        return Err(exc);
    }

    let mut log = String::new();
    if shared_memory.unload_last {
        if in_use > 0 {
            log = format!(
                "{} other job(s) are attached to the shared memory segment, will not remove it.\n",
                in_use
            );
        } else {
            log = "No other jobs are attached to the shared memory segment, removing it.\n"
                .to_string();
            if let Err(exc) = sharedmemory_l237_sharedmemory_clean(shared_memory) {
                let _ = sharedmemory_l237_sharedmemory_clean(shared_memory);
                return Err(exc);
            }
        }
    }
    Ok(log)
}

#[doc = "Original `SharedMemory::Allocate` at STAR/source/SharedMemory.cpp:68. Args: shmSize: size_t"]
pub fn sharedmemory_l68_sharedmemory_allocate(
    shared_memory: &mut crate::shared_memory::SharedMemory,
    shm_size: usize,
) -> Result<(), crate::shared_memory::SharedMemoryException> {
    shared_memory.exception = crate::shared_memory::SharedMemoryException::default();

    if !shared_memory.needs_allocation {
        shared_memory.exception = crate::shared_memory::SharedMemoryException {
            error_code: 3,
            error_detail: 0,
        };
        return Err(shared_memory.exception.clone());
    }

    sharedmemory_l102_sharedmemory_createandinitsharedobject(shared_memory, shm_size)?;
    if shared_memory.exception.error_code != 0 && shared_memory.exception.error_code != SHM_EEXISTS
    {
        return Err(shared_memory.exception.clone());
    }

    shared_memory.exception = crate::shared_memory::SharedMemoryException::default();
    sharedmemory_l134_sharedmemory_openifexists(shared_memory)?;
    shared_memory.is_allocator = true;
    Ok(())
}

#[doc = "Original `SharedMemory::GetPosixObjectKey` at STAR/source/SharedMemory.cpp:87. Args: "]
pub fn sharedmemory_l87_sharedmemory_getposixobjectkey(
    shared_memory: &crate::shared_memory::SharedMemory,
) -> String {
    format!("/{}", shared_memory.key)
}

#[doc = "Original `SharedMemory::CounterName` at STAR/source/SharedMemory.cpp:94. Args: "]
pub fn sharedmemory_l94_sharedmemory_countername(
    shared_memory: &crate::shared_memory::SharedMemory,
) -> String {
    format!("/shared_use_counter{}", shared_memory.key)
}

#[doc = "Original `SharedMemory::CreateAndInitSharedObject` at STAR/source/SharedMemory.cpp:102. Args: shmSize: size_t"]
pub fn sharedmemory_l102_sharedmemory_createandinitsharedobject(
    shared_memory: &mut crate::shared_memory::SharedMemory,
    shm_size: usize,
) -> Result<(), crate::shared_memory::SharedMemoryException> {
    if shared_memory.shm_id != -1 {
        shared_memory.exception = crate::shared_memory::SharedMemoryException {
            error_code: SHM_EEXISTS,
            error_detail: 0,
        };
        return Ok(());
    }
    shared_memory.shm_id = shared_memory.key;
    shared_memory.length = shm_size + std::mem::size_of::<usize>();
    Ok(())
}

#[doc = "Original `SharedMemory::OpenIfExists` at STAR/source/SharedMemory.cpp:134. Args: "]
pub fn sharedmemory_l134_sharedmemory_openifexists(
    shared_memory: &mut crate::shared_memory::SharedMemory,
) -> Result<bool, crate::shared_memory::SharedMemoryException> {
    let exists = shared_memory.shm_id >= 0;
    if exists {
        sharedmemory_l168_sharedmemory_mapsharedobjecttomemory(shared_memory)?;
        shared_memory.needs_allocation = false;
    }
    Ok(exists)
}

#[doc = "Original `SharedMemory::GetSharedObjectInfo` at STAR/source/SharedMemory.cpp:157. Args: "]
pub fn sharedmemory_l157_sharedmemory_getsharedobjectinfo(
    shared_memory: &crate::shared_memory::SharedMemory,
) -> Result<usize, crate::shared_memory::SharedMemoryException> {
    if shared_memory.shm_id == -1 {
        Err(crate::shared_memory::SharedMemoryException {
            error_code: SHM_EOPENFAILED,
            error_detail: 0,
        })
    } else {
        Ok(shared_memory.length)
    }
}

#[doc = "Original `SharedMemory::MapSharedObjectToMemory` at STAR/source/SharedMemory.cpp:168. Args: "]
pub fn sharedmemory_l168_sharedmemory_mapsharedobjecttomemory(
    shared_memory: &mut crate::shared_memory::SharedMemory,
) -> Result<(), crate::shared_memory::SharedMemoryException> {
    if shared_memory.shm_id == -1 {
        shared_memory.exception = crate::shared_memory::SharedMemoryException {
            error_code: SHM_EMAPFAILED,
            error_detail: 0,
        };
        return Err(shared_memory.exception.clone());
    }
    shared_memory.mapped = true;
    Ok(())
}

#[doc = "Original `SharedMemory::Close` at STAR/source/SharedMemory.cpp:191. Args: "]
pub fn sharedmemory_l191_sharedmemory_close(
    shared_memory: &mut crate::shared_memory::SharedMemory,
) -> Result<(), crate::shared_memory::SharedMemoryException> {
    shared_memory.mapped = false;
    shared_memory.shm_id = -1;
    Ok(())
}

#[doc = "Original `SharedMemory::Unlink` at STAR/source/SharedMemory.cpp:219. Args: "]
pub fn sharedmemory_l219_sharedmemory_unlink(
    shared_memory: &mut crate::shared_memory::SharedMemory,
) -> Result<bool, crate::shared_memory::SharedMemoryException> {
    if !shared_memory.needs_allocation {
        shared_memory.needs_allocation = true;
        Ok(true)
    } else {
        Ok(false)
    }
}

#[doc = "Original `SharedMemory::Clean` at STAR/source/SharedMemory.cpp:237. Args: "]
pub fn sharedmemory_l237_sharedmemory_clean(
    shared_memory: &mut crate::shared_memory::SharedMemory,
) -> Result<(), crate::shared_memory::SharedMemoryException> {
    sharedmemory_l191_sharedmemory_close(shared_memory)?;
    sharedmemory_l219_sharedmemory_unlink(shared_memory)?;
    sharedmemory_l269_sharedmemory_removesharedcounter(shared_memory)?;
    shared_memory.clean_count += 1;
    Ok(())
}

#[doc = "Original `SharedMemory::EnsureCounter` at STAR/source/SharedMemory.cpp:244. Args: "]
pub fn sharedmemory_l244_sharedmemory_ensurecounter(
    shared_memory: &mut crate::shared_memory::SharedMemory,
) -> Result<(), crate::shared_memory::SharedMemoryException> {
    if shared_memory.shared_counter_id < 0 {
        shared_memory.shared_counter_id = shared_memory.counter_key;
    }
    if !shared_memory.counter_mem_attached {
        shared_memory.counter_mem_attached = true;
    }
    Ok(())
}

#[doc = "Original `SharedMemory::RemoveSharedCounter` at STAR/source/SharedMemory.cpp:269. Args: "]
pub fn sharedmemory_l269_sharedmemory_removesharedcounter(
    shared_memory: &mut crate::shared_memory::SharedMemory,
) -> Result<(), crate::shared_memory::SharedMemoryException> {
    if shared_memory.shared_counter_id == -1 {
        shared_memory.exception = crate::shared_memory::SharedMemoryException {
            error_code: SHM_ECOUNTERREMOVE,
            error_detail: 0,
        };
        return Err(shared_memory.exception.clone());
    }
    shared_memory.shared_counter_id = -1;
    shared_memory.counter_mem_attached = false;
    Ok(())
}

#[doc = "Original `SharedMemory::SharedObjectsUseCount` at STAR/source/SharedMemory.cpp:277. Args: "]
pub fn sharedmemory_l277_sharedmemory_sharedobjectsusecount(
    shared_memory: &mut crate::shared_memory::SharedMemory,
) -> Result<i32, crate::shared_memory::SharedMemoryException> {
    sharedmemory_l244_sharedmemory_ensurecounter(shared_memory)?;
    if shared_memory.shared_counter_id != -1 {
        Ok(shared_memory.shared_objects_use_count_value)
    } else {
        Ok(-1)
    }
}
